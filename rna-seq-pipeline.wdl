# ENCODE DCC RNA-seq pipeline
# Maintainer: Otto Jolanki
# Modified for 4DN by Soo Lee

#CAPER docker quay.io/encode-dcc/rna-seq-pipeline:v1.1
#CAPER singularity docker://quay.io/encode-dcc/rna-seq-pipeline:v1.1
#CROO out_def https://storage.googleapis.com/encode-pipeline-output-definition/bulkrna.output_definition.json

workflow rna {
    # endedness: paired or single
    String endedness
    # fastqs_R1: fastq.gz files for Read1 (only these if single-ended)
    Array[Array[File]] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[Array[File]] fastqs_R2 = []
    # bamroot: root name for output bams. For example foo_bar will
    # create foo_bar_genome.bam and foo_bar_anno.bam
    String bamroot = "out"
    # strandedness: is the library strand specific (stranded or unstranded)
    String strandedness
    # strandedness_direction (forward, reverse, unstranded)
    String strandedness_direction
    # chrom_sizes: chromosome sizes file
    File chrom_sizes

    Int cpu = 16
    Int memGB = 64

    ## task level variables that are defined globally to make them visible to DNANexus UI

    # ALIGN
    # index: aligner index archive (tar.gz)
    File align_index
    Int align_ncpus = cpu
    Int align_ramGB = memGB
    String? align_disk

    # BAM_TO_SIGNALS

    Int bam_to_signals_ncpus = cpu
    Int bam_to_signals_ramGB = memGB
    String? bam_to_signals_disk

    # RSEM_QUANT

    # rsem_index: RSEM index archive (tar.gz)
    File rsem_index
    # rnd_seed: random seed used for rsem
    Int rnd_seed = 12345
    Int rsem_ncpus = cpu
    Int rsem_ramGB = memGB
    String? rsem_disk

    # RNA_QC

    File rna_qc_tr_id_to_gene_type_tsv
    String? rna_qc_disk

    # MAD_QC

    String? mad_qc_disk

    ## WORKFLOW BEGINS

    # dummy variable value for the single-ended case
    Array[Array[File]] fastqs_R2_ = if (endedness == "single") then fastqs_R1 else fastqs_R2

    scatter (i in range(length(fastqs_R1))) {
        call align { input:
            endedness = endedness,
            fastqs_R1 = fastqs_R1[i],
            fastqs_R2 = fastqs_R2_[i],
            index = align_index,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = align_disk,
        }

        call bam_to_signals { input:
            input_bam = align.genomebam,
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "rep"+(i+1)+bamroot+"_genome",
            ncpus = bam_to_signals_ncpus,
            ramGB = bam_to_signals_ramGB,
            disks = bam_to_signals_disk,
        }

        call rsem_quant { input:
            rsem_index = rsem_index,
            rnd_seed = rnd_seed,
            anno_bam = align.annobam,
            endedness = endedness,
            read_strand = strandedness_direction,
            ncpus = rsem_ncpus,
            ramGB = rsem_ramGB,
            disks = rsem_disk,
        }
    }

    # if there are exactly two replicates, calculate the madQC metrics and draw a plot
    if (length(fastqs_R1) == 2) {
        call mad_qc { input:
            quants1 = rsem_quant.genes_results[0],
            quants2 = rsem_quant.genes_results[1],
            disks = mad_qc_disk,
        }
    }

    scatter (i in range(length(align.annobam))) {
        call rna_qc { input:
            input_bam = align.annobam[i],
            tr_id_to_gene_type_tsv = rna_qc_tr_id_to_gene_type_tsv,
            output_filename = "rep"+(i+1)+bamroot+"_qc.json",
            disks = rna_qc_disk,
        }
    }


    if(strandedness == 'unstranded'){
        File? unstranded_bw = bam_to_signals.unique[0][0]
    }
    if(strandedness != 'unstranded'){
        File? minus_bw = bam_to_signals.unique[0][0]
        File? plus_bw = bam_to_signals.unique[0][1]
    }
    
    output {
        File outbam = align.genomebam[0]
        File outbw = select_first([unstranded_bw, ''])
        File minusbw = select_first([minus_bw, ''])
        File plusbw = select_first([plus_bw, ''])
        File gene_expression = rsem_quant.genes_results[0]
        File isoform_expression = rsem_quant.isoforms_results[0]
        File genome_flagstat_qc_json = align.genome_flagstat_json[0]
        File anno_flagstat_qc_json = align.anno_flagstat_json[0]
        File align_log_json = align.log_json[0]
        File qc_json = rna_qc.rnaQC[0]
    }
}

## tasks
task align {
    Array[File] fastqs_R1
    Array[File] fastqs_R2
    String endedness
    File index
    String bamroot
    Int ncpus
    Int ramGB
    String? disks

    command {
        python3 $(which align.py) \
            --fastqs_R1 ${sep=' ' fastqs_R1} \
            --fastqs_R2 ${sep=' ' fastqs_R2} \
            --endedness ${endedness} \
            --index ${index} \
            ${"--bamroot " + bamroot} \
            ${"--ncpus " + ncpus} \
            ${"--ramGB " + ramGB}
    }

    output {
        File genomebam = glob("*_genome.bam")[0]
        File annobam = glob("*_anno.bam")[0]
        File genome_flagstat = glob("*_genome_flagstat.txt")[0]
        File anno_flagstat = glob("*_anno_flagstat.txt")[0]
        File log = glob("*_Log.final.out")[0]
        File genome_flagstat_json = glob("*_genome_flagstat.json")[0]
        File anno_flagstat_json = glob("*_anno_flagstat.json")[0]
        File log_json = glob("*_Log.final.json")[0]
        File python_log = glob("align.log")[0]
    }

    runtime {
      cpu: ncpus
      memory: "${ramGB} GB"
      disks : select_first([disks,"local-disk 100 SSD"])
      docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}

task  bam_to_signals {
    File input_bam
    File chrom_sizes
    String strandedness
    String bamroot
    Int ncpus
    Int ramGB
    String? disks


    command {
        python3 $(which bam_to_signals.py) \
            --bamfile ${input_bam} \
            --chrom_sizes ${chrom_sizes} \
            --strandedness ${strandedness} \
            --bamroot ${bamroot}
    }

    output {
        Array[File] unique = glob("*niq.bw")
        Array[File] all = glob("*ll.bw")
        File python_log = glob("bam_to_signals.log")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
        docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}

task rsem_quant {
    File rsem_index
    File anno_bam
    String endedness
    String read_strand
    Int rnd_seed
    Int ncpus
    Int ramGB
    String? disks

    command {
        python3 $(which rsem_quant.py) \
            --rsem_index ${rsem_index} \
            --anno_bam ${anno_bam} \
            --endedness ${endedness} \
            --read_strand ${read_strand} \
            --rnd_seed ${rnd_seed} \
            --ncpus ${ncpus} \
            --ramGB ${ramGB}
    }

    output {
        File genes_results = glob("*.genes.results")[0]
        File isoforms_results = glob("*.isoforms.results")[0]
        File python_log = glob("rsem_quant.log")[0]
        File number_of_genes = glob("*_number_of_genes_detected.json")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
        docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}

task kallisto {
    Array[File] fastqs_R1
    Array[File] fastqs_R2
    File kallisto_index
    String endedness
    String strandedness_direction
    Int number_of_threads
    Int ramGB
    String out_prefix
    Int? fragment_length
    Float? sd_of_fragment_length
    String? disks

    command {
        python3 $(which kallisto_quant.py) \
            --fastqs_R1 ${sep=' ' fastqs_R1} \
            --fastqs_R2 ${sep=' ' fastqs_R2} \
            --number_of_threads ${number_of_threads} \
            --strandedness ${strandedness_direction} \
            --path_to_index ${kallisto_index} \
            --endedness ${endedness} \
            ${"--fragment_length " + fragment_length} \
            ${"--sd_of_fragment_length " + sd_of_fragment_length} \
            ${"--out_prefix " + out_prefix}
    }

    output {
        File quants = glob("kallisto_out/*_abundance.tsv")[0]
        File python_log = glob("kallisto_quant.log")[0]
    }

    runtime {
        cpu: number_of_threads
        memory: "${ramGB} GB"
        disks: select_first([disks, "local-disk 100 SSD"])
        docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}

task mad_qc {
    File quants1
    File quants2
    String? disks

    command {
        python3 $(which mad_qc.py) \
            --quants1 ${quants1} \
            --quants2 ${quants2} \
            --MAD_R_path $(which MAD.R)
    }

    output {
        File madQCplot = glob("*_mad_plot.png")[0]
        File madQCmetrics = glob("*_mad_qc_metrics.json")[0]
        File python_log = glob("mad_qc.log")[0]
    }

    runtime {
        cpu: 2
        memory: "3400 MB"
        disks: select_first([disks,"local-disk 100 SSD"])
        docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}

task rna_qc {
    File input_bam
    File tr_id_to_gene_type_tsv
    String output_filename
    String? disks

    command {
        python3 $(which rna_qc.py) \
            --input_bam ${input_bam} \
            --tr_id_to_gene_type_tsv ${tr_id_to_gene_type_tsv} \
            --output_filename ${output_filename}
    }

    output {
        File rnaQC = glob("*_qc.json")[0]
        File python_log = glob("rna_qc.log")[0]
    }

    runtime {
        cpu: 2
        memory: "1024 MB"
        disks: select_first([disks, "local-disk 100 SSD"])
        docker: "4dndcic/encode-rnaseq:v1.1_dcic"
    }
}
