nextflow.enable.dsl = 2

workflow {
    input_data = Channel.fromPath(params.samplesheet)
        | splitCsv(sep: '\t', header: true)
        | map { row ->
            [
                row.ID_family,
                [file(row.BAM_c), file(row.BAM_c_i)],  // child BAM + index
                [file(row.BAM_m), file(row.BAM_m_i)],  // mother BAM + index
                [file(row.BAM_f), file(row.BAM_f_i)],  // father BAM + index
                [file(row.VCF_m), file(row.VCF_m_tbi)],  // mother VCF + index
                [file(row.VCF_f), file(row.VCF_f_tbi)],  // father VCF + index
                [file(row.VCF_c), file(row.VCF_c_tbi)]   // child VCF + index
            ]
        }
        | map { id, bam_c, bam_m, bam_f, vcf_m, vcf_f, vcf_c -> 
            tuple(
                id,
                bam_c[0], bam_c[1],
                bam_m[0], bam_m[1],
                bam_f[0], bam_f[1],    
                vcf_m[0], vcf_m[1],
                vcf_f[0], vcf_f[1],
                vcf_c[0], vcf_c[1]
            )
        }

    input_data.view { it -> "Списки образцов на вход: ${it}" }

    merged_vcf_input = input_data.map { tuple(it[0], it[7], it[8], it[9], it[10], it[11], it[12]) }

    merged_vcf_ch = MergeVCF(merged_vcf_input)

    reference_ch = Channel.fromPath(params.reference)
    ped_ch = Channel.fromPath(params.ped)

    PhaseVCF_input = merged_vcf_ch
        .combine(reference_ch)
        .combine(ped_ch)
        .join(
            input_data.map { id, bams_c0, bams_c1, bams_m0, bams_m1, bams_f0, bams_f1, vcf_m0, vcf_m1, vcf_f0, vcf_f1, vcf_c0, vcf_c1 -> 
                tuple(
                    id, 
                    bams_c0, bams_c1,  // child.bam, child.bai
                    bams_m0, bams_m1,  // mother.bam, mother.bai
                    bams_f0, bams_f1   // father.bam, father.bai
                ) 
            }, 
            by: 0
        )
        .map { id, merged_vcf, ref, ped, c_bam, c_bai, m_bam, m_bai, f_bam, f_bai -> 
            tuple(
                id, 
                merged_vcf,
                ref, ped,
                c_bam, c_bai,
                m_bam, m_bai,
                f_bam, f_bai
            ) 
        }

    PhaseVCF(PhaseVCF_input, "/media/HEAP-EPI/etcetera_mi/References/hg38_btk/GRCh38.d1.vd1.fa.fai")
    PhaseVCF_input.view { it -> "PhaseVCF списки на вход: ${it}" }
}

process MergeVCF {
    tag "${id}"
    publishDir "${params.outdir}/merged_vcf/${id}", mode: 'copy'
    container "${ workflow.containerEngine == 'docker' ? 'eod-tools.med-gen.ru/icr_pipeline:latest' : 
    'eod-tools.med-gen.ru/icr_pipeline:latest' }"

    input:
    tuple val(id),
      path(vcf_m, stageAs: "mother.vcf.gz"),
      path(vcf_m_tbi, stageAs: "mother.vcf.gz.tbi"),
      path(vcf_f, stageAs: "father.vcf.gz"),
      path(vcf_f_tbi, stageAs: "father.vcf.gz.tbi"),
      path(vcf_c, stageAs: "child.vcf.gz"),
      path(vcf_c_tbi, stageAs: "child.vcf.gz.tbi")


    output:
    tuple val(id), path("${id}_merged.vcf")

    script:
    """
    bcftools merge -o "${id}_merged.vcf" mother.vcf.gz father.vcf.gz child.vcf.gz
    """
}

process PhaseVCF {
    tag "${id}"
    publishDir "${params.outdir}/phased_vcf/${id}", mode: 'copy'
    container "${ workflow.containerEngine == 'docker' ? 'hkubal/clair3-nova:latest' : 
    'hkubal/clair3-nova:latest' }"

    input:
    tuple val(id), 
          path(merged_vcf),    // Объединенный VCF
          path(reference),     // Референс
          path(ped),           // PED 
          path(c_bam),         // Child BAM
          path(c_bai),         // Child BAI
          path(m_bam),         // Mother BAM
          path(m_bai),         // Mother BAI
          path(f_bam),         // Father BAM
          path(f_bai)          // Father BAI
    path fai_index  

    output:
    tuple val(id), 
          path("${id}.phased.vcf"), 
          path("${id}_stats.tsv"), 
          path("${id}_phased.gtf")

    script:
    """
    export FAI_index=${fai_index}
    ${params.whatshap} phase -o ${id}.phased.vcf --ped ${ped} --reference ${reference} --tag HP ${merged_vcf} ${c_bam} ${m_bam} ${f_bam}
    ${params.whatshap} stats --gtf=${id}_phased.gtf --tsv=${id}_stats.tsv ${id}.phased.vcf
    """
}
