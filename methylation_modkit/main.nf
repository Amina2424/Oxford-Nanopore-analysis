nextflow.enable.dsl = 2

process modkit_pileup {

    container "eod-tools.med-gen.ru/icr_pipeline:latest"

    publishDir "${params.results_dir}/pileups", mode: 'copy'

    tag "${id}"
    
    input:
    tuple val(id), path(bam), path(bai)
    path reference
    path fai

    output:
    tuple val(id), path("${id}_pileup.bed")
    
    script:
    """
    modkit pileup ${bam} \
        --cpg \
        --ignore h \
        ${id}_pileup.bed \
	--ref ${reference} \
        --combine-strands \
        -t 35
    """
}

process bgzip_tabix {

    container "${ workflow.containerEngine == 'docker' ? 'eod-tools.med-gen.ru/icr_pipeline:latest' : 
    'eod-tools.med-gen.ru/icr_pipeline:latest' }"
    publishDir "${params.results_dir}/compressed", mode: 'copy'

    tag "${id}"
    
    input:
    tuple val(id), path(bed)
    
    output:
    tuple val(id), path("${id}_pileup.bed.gz"), path("${id}_pileup.bed.gz.tbi")
    
    script:
    """
    bgzip -c ${id}_pileup.bed > ${id}_pileup.bed.gz
    tabix -p bed ${id}_pileup.bed.gz
    """
}

process merge_all {

    container "${ workflow.containerEngine == 'docker' ? 'eod-tools.med-gen.ru/icr_pipeline:latest' : 
    'eod-tools.med-gen.ru/icr_pipeline:latest' }"
    publishDir "${params.results_dir}/merged", mode: 'copy'

    input:
    path bedgzs
    path genome_size
    path tbis

    output:
    path "merged_bedmethyl.bed.gz"
    
    script:
    """

    modkit bedmethyl merge ${bedgzs} --out-bed merged_bedmethyl.bed --genome-sizes ${genome_size} --threads 35

    bgzip merged_bedmethyl.bed
    tabix merged_bedmethyl.bed.gz

    """
}

process roi_stats {

    container "${ workflow.containerEngine == 'docker' ? 'eod-tools.med-gen.ru/icr_pipeline:latest' : 
    'eod-tools.med-gen.ru/icr_pipeline:latest' }"
    publishDir "${params.results_dir}/stats", mode: 'copy'
   
    input:
    path merged_bed
    path roi_bed
    
    output:
    path "ctrl_gDMRs_stats.tsv"
    
    script:
    """
    
    tabix ${merged_bed}
    modkit stats \
        -o ctrl_gDMRs_stats.tsv \
        --regions ${roi_bed} \
        ${merged_bed}
    """
}

workflow {
    ch_reference = "/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa"
    ch_genome_sizes = Channel.fromPath(params.genome_sizes)
    ch_reference_fai = "/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa.fai"

    samples_ch = Channel.fromPath(params.samplesheet)
        | splitCsv( sep: "\t", header: true, strip: true)
        | map { row -> tuple(
            row.ID,
            //file(row.BAM, checkIfExists: true*/),
            //file(row.BAI, checkIfExists: true)
	    row.BAM,
	    row.BAI
        ) 
        }
//    samples_ch.view()

    ch_pileup = modkit_pileup(samples_ch, ch_reference, ch_reference_fai)
/*        .set { pileup_ch }

    bgzip = bgzip_tabix(pileup_ch)
        
    gz = bgzip.map { id, gz, tbi -> gz } 
        .collect()

    tbi = bgzip.map { id, gz, tbi -> tbi }
        .collect()
 
    merge_all(gz,"/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa.fai", tbi)

    roi_stats(merge_all.out, params.roi_bed)*/

}
