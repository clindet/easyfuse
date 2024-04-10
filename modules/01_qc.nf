
process FASTP {
    cpus 6
    memory "8g"
    tag "${name}"

    conda '/public/home/lijf/env/miniconda3/envs/easyfuse-qc'

    input:
      tuple val(name), path(fastq1), path(fastq2)

    output:
      tuple val("${name}"), path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), emit: trimmed_fastq

    """
    fastp \
        -i ${fastq1} \
        -I ${fastq2} \
        -o trimmed_R1.fastq.gz \
        -O trimmed_R2.fastq.gz \
        --thread ${task.cpus}
    """
}
