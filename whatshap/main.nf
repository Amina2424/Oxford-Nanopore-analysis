nextflow.enable.dsl = 2

Channel.fromPath(params.samplesheet)
    .splitCsv(sep: "\t", header: true)
    .map { row -> 
        tuple(
            row.ID,
            file(row.BAM, checkExists: true),
            file(row.BAI, checkExists: true),
            file(row.VCF, checkExists: true),
            file(row.VCF_tbi, checkExists: true)
        ) 
    }
    .set { samples_ch }

workflow {

    ch_reference = "/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa"
    ch_reference_fai = "/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa.fai"
    bed = "/media/HEAP-EPI/etcetera_mi/hg_38_res/candidate_dmrs/reference/merged_candidates.bed"

    WhatshapHaplotag(samples_ch, ch_reference, ch_reference_fai)
     .set { haplotagged_ch }
    
    modkit_pileup(haplotagged_ch, ch_reference, ch_reference_fai, bed)
     .set { pileup_ch }
}

process WhatshapHaplotag {
    tag "${id}"
    publishDir "${params.outdir}/haplotagged_bams/${id}", mode: 'copy'
    container "eod-tools.med-gen.ru/icr_pipeline:latest"

    input:
    tuple val(id), 
          path(bam),
          path(bai),
          path(vcf),
          path(vcf_tbi)
    path reference
    path fai

    output:
    tuple val(id), 
          path("hp_${id}.bam"), 
          path("hp_${id}.bam.bai")
    
    script:
    """
    export TBI_index=${vcf_tbi}
    ${params.whatshap} haplotag \\
        -o hp_${id}.bam \\
        ${vcf} \\
        ${bam} \\
        --reference ${reference} \\
        --output-threads 35

    samtools index hp_${id}.bam
    """
}

process modkit_pileup {
    tag "${id}"
    publishDir "${params.outdir}/pileups/${id}", mode: 'copy'
    container "eod-tools.med-gen.ru/icr_pipeline:latest"

    input:
    tuple val(id), 
          path(bam),
          path(bai)
    path reference
    path fai
    path roi_bed  // Add this input for your regions of interest BED file

    output:
    tuple val(id), path("${id}_temp/*")
    tuple val(id), path("${id}_stats/*")

    script:
    """
    mkdir -p ${id}_temp
    mkdir -p ${id}_stats
    
   
    export FAI_index=${fai}
    ${params.modkit} pileup ${bam} \\
        --cpg \\
        --ref ${reference} \\
        --ignore h \\
        -t 35 \\
        --partition-tag HP \\
        --combine-strands \\
        --filter-threshold C:0.7 \\
        --mod-threshold m:0.9 \\
        ${id}_temp/
    
    for bed_file in ${id}_temp/*.bed; do
        if [ -s "\$bed_file" ]; then

            bgzip -f -@ 10 "\$bed_file"
            
            tabix -f -p bed "\${bed_file}.gz"
            
            base_name=\$(basename "\$bed_file" .bed)
            ${params.modkit} stats --regions ${roi_bed} "\${bed_file}.gz" -o "${id}_stats/\${base_name}_stats.csv"
        fi
    done
    """
}
