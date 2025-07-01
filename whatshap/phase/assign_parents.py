#!/usr/bin/env python3

from cyvcf2 import VCF

# Загрузка VCF-файлов
child_vcf = VCF("/media/BTK/2024/05/5-ONT/ont2/000007030520-ONT/hg38/vcf/Clair3_phased.vcf.gz")
mother_vcf = VCF("/media/BTK/2024/05/5-ONT/ont2/000007030530-ONT/hg38/vcf/Clair3_phased.vcf.gz")
father_vcf = VCF("/media/BTK/2024/05/5-ONT/ont2/000007030540-ONT/hg38/vcf/Clair3_phased.vcf.gz")

# Определение имён сэмплов
mother_id = mother_vcf.samples[0]
father_id = father_vcf.samples[0]
child_id = child_vcf.samples[0]

# Индексация вариантов по позиции
def index_variants(vcf_reader, sample_name):
    index = {}
    sample_idx = vcf_reader.samples.index(sample_name)
    for rec in vcf_reader:
        if len(rec.ALT) != 1:
            continue  # Только биалельные варианты
        gt = rec.genotypes[sample_idx][:2]
        if None in gt:
            continue
        key = (rec.CHROM, rec.POS)
        index[key] = {
            "ref": rec.REF,
            "alt": rec.ALT[0],
            "gt": gt
        }
    return index

# Индексируем материнские и отцовские варианты
mother_variants = index_variants(mother_vcf, mother_id)
father_variants = index_variants(father_vcf, father_id)

# Заголовок
print("CHROM\tPOS\tChild_GT\tPhase0\tPhase1\tHP1\tHP2")

# Обход вариантов в ребёнке
child_sample_idx = child_vcf.samples.index(child_id)

for rec in child_vcf:
    if len(rec.ALT) != 1:
        continue

    if not rec.gt_phases[child_sample_idx]:  # Проверка фазировки
        continue

    gt = rec.genotypes[child_sample_idx][:2]
    if gt not in [(0, 1), (1, 0)]:
        continue

    key = (rec.CHROM, rec.POS)
    if key not in mother_variants or key not in father_variants:
        continue

    ref = rec.REF
    alt = rec.ALT[0]
    alleles = [ref, alt]

    phase0 = alleles[gt[0]]
    phase1 = alleles[gt[1]]

    mother_gt = mother_variants[key]["gt"]
    father_gt = father_variants[key]["gt"]
    mother_alleles = set([alleles[x] for x in mother_gt])
    father_alleles = set([alleles[x] for x in father_gt])

    hp1 = hp2 = "?"

    if phase0 in mother_alleles and phase1 in father_alleles:
        hp1 = "maternal"
        hp2 = "paternal"
    elif phase1 in mother_alleles and phase0 in father_alleles:
        hp1 = "paternal"
        hp2 = "maternal"

    print(f"{rec.CHROM}\t{rec.POS}\t{gt[0]}|{gt[1]}\t{phase0}\t{phase1}\t{hp1}\t{hp2}")
