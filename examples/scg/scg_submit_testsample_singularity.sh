#!/bin/bash

# do not touch these settings
# number of tasks and nodes are fixed at 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1

# job name for pipeline
# this name will appear when you monitor jobs with "squeue -u $USER"
#SBATCH --job-name=RNA-SEQ-testsample
# walltime for your job
# give long time enough to finish your pipeline
#SBATCH --time=12:00:00

# total amount of memory
# depends on the size of your FASTQs
# do not request too much memory
# cluster will not accept your job
#SBATCH --mem=20G

# max number of cpus for each pipeline
# should be >= NUM_CONCURRENT_TASK x "rna.align_ncpus"
#            + NUM_CONCURRENT_TASK x "rna.kallisto_number_of_threads"
# in input JSON file. This is the worst case maximum number
# of cpus the pipeline may use. 
#SBATCH --cpus-per-task=8

# email notification for job status
#SBATCH --mail-type=END,FAIL

# load java module if it exists
module load java || true

# use input JSON for a small test sample
# make an input JSON for your own sample
# start from any of two templates for single-ended and paired-ended samples
# (examples/template_se.json, examples/template_pe.json)
INPUT=examples/scg/sct_testsample_input.json

# If this pipeline fails, then use this metadata JSON file to resume a failed pipeline from where it left 
# See details in /utils/resumer/README.md
PIPELINE_METADATA=metadata.json

# limit number of concurrent tasks
#  we recommend to use a number of replicates here
#  so that all replicates are processed in parellel at the same time.
#  make sure that resource settings in your input JSON file
#  are consistent with SBATCH resource settings (--mem, --cpus-per-task) 
#  in this script
NUM_CONCURRENT_TASK=2

# run pipeline
#  you can monitor your jobs with "squeue -u $USER"
java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=singularity \
-Dbackend.providers.singularity.config.concurrent-job-limit=${NUM_CONCURRENT_TASK} \
$HOME/cromwell-35.jar run rna-seq-pipeline.wdl -i ${INPUT} -o workflow_opts/scg.json -m ${PIPELINE_METADATA}