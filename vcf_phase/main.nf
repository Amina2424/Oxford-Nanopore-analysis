nextflow.enable.dsl = 2

workflow {
    input_data = Channel.fromPath(params.samplesheet)
        | splitCsv(sep: '\t', header: true)
        | map { row -> tuple(row.ID, file(row.VCF), file(row.VCF_tbi)) }
    
    reference_ch = Channel.fromPath(params.reference)
    fai_ch = Channel.fromPath("/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa.fai")

    PhaseVCF_input = input_data
        .combine(reference_ch)
        .combine(fai_ch)
    
    PhaseVCF_input.view { id, vcf, tbi, ref, fai -> 
        "PhaseVCF входные данные: $id, VCF: ${vcf.name}, Reference: ${ref.name}, FAI: ${fai.name}" 
    }

    PhaseVCF(PhaseVCF_input)
}

process PhaseVCF {
    tag "${id}"
    publishDir "${params.outdir}/phased_vcf/${id}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".gtf")) "gtf/$filename"
            else filename
        }
    
    container "eod-tools.med-gen.ru/icr_pipeline:latest"

    input:
    tuple val(id), 
          path(vcf),
          path(vcf_tbi),
          path(reference),
          path(fai_index)

    output:
    tuple val(id), 
          path("${id}_stats.tsv"), 
          path("${id}_phased.gtf")

    script:
    """
    export FAI_index=${fai_index}
    ${params.whatshap} stats \
        --gtf=${id}_phased.gtf \
        --tsv=${id}_stats.tsv \
        ${vcf}
    """
}
