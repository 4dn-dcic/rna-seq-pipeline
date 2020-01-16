workflow mad_qc {
	Array[File] quantfiles
	String? mad_qc_disk
    	Array[String]? sampleids
	Array[String]? urls

	call mqc { input:
		quants = quantfiles,
		disks = mad_qc_disk,
		ids = sampleids,
		links = urls
	}
}

task mqc {
	Array[File] quants
	String? disks
	Array[String]? ids
	Array[String]? links

	command {
	python3 $(which mad_qc.py) \
			--quants ${sep=' ' quants} \
			--MAD_R_path $(which MAD.R) \
			--sample_ids ${sep=' ' ids} \
			--quant_urls ${sep=' ' links}
	}

	output {
		File report_zip = glob("mad_qc_report.zip")[0]
		File madQCmetrics = glob("*_mad_qc_metrics.json")[0]
	}

	runtime {
		cpu: 2
		memory: "3400 MB"
		disks: select_first([disks,"local-disk 100 SSD"])
		docker: "4dndcic/encode-rnaseq:v1.1_dcic_2"

	}
}
