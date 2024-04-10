
process FUSION_CATCHER_INDEX {
    cpus 10
    memory "32g"
    tag "${name}"

    //conda ("${baseDir}/environments/fusioncatcher.yml")
    conda "/public/home/lijf/env/miniconda3/envs/easyfuse-fusioncatcher"

    script:
    """
    download-human-db.sh
    """
}


process FUSION_CATCHER {
    cpus 10
    memory "32g"
    tag "${name}"

    //conda ("${baseDir}/environments/fusioncatcher.yml")
    conda "/public/home/lijf/env/miniconda3/envs/easyfuse-fusioncatcher"
    publishDir "${params.output}/${name}/fusioncatcher", mode: 'copy'

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("summary_candidate_fusions.txt"), path("final-list_candidate-fusion-genes.txt"), emit: fusions

    script:
    """
    
    fusioncatcher \
        --data ${params.fusioncatcher_index} \
        --input ${fastq1},${fastq2} \
        --output . \
        -p ${task.cpus}
    """
}

process STAR_FUSION {
    cpus 10
    memory "32g"
    tag "${name}"

    //conda ("${baseDir}/environments/starfusion.yml")
    conda "/public/home/lijf/env/miniconda3/envs/easyfuse-starfusion"
    publishDir "${params.output}/${name}/star_fusion", mode: 'copy'

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("star-fusion.fusion_predictions.tsv"), emit: fusions

    script:
    """
    STAR-Fusion \
        --left_fq ${fastq1} \
        --right_fq ${fastq2} \
        --genome_lib_dir ${params.starfusion_index} \
        --CPU ${task.cpus} \
        --output_dir .
    """
}

process ARRIBA {
    cpus 10
    memory "20g"
    tag "${name}"

    //conda ("${baseDir}/environments/arriba.yml")
    conda "/public/home/lijf/env/miniconda3/envs/easyfuse-arriba"
    publishDir "${params.output}/${name}/arriba", mode: 'copy'

    input:
      tuple val(name), path(bam)

    output:
      tuple val("${name}"), path("fusions.tsv"), path("fusions.discarded.tsv"), emit: fusions

    script:
    """
    arriba \
        -x ${bam} \
        -g ${params.gtf} \
        -a ${params.fasta} \
        -o fusions.tsv \
        -O fusions.discarded.tsv \
        -f blacklist
    """
}
