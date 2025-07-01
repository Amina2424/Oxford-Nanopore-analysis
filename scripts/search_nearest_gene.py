#!/usr/bin/env python3

import pyensembl
import pandas as pd

ensembl = pyensembl.EnsemblRelease(release=110)  # релиз
print("Загружены аннотации для сборки:", pyensembl.ensembl_grch38)

# Функция для поиска ближайшего транскрипта
def find_nearest_transcript(chrom, start, end):
    chrom = chrom.replace("chr", "")  # у энсэмбл без хром
    transcripts = ensembl.transcripts_at_locus(chrom, start, end)
    if not transcripts:
        return None
    # Возвращаем первый транскрипт и его цепь
    return transcripts[0], transcripts[0].strand, transcripts[0].biotype

file_path = "merged_candidates.bed"
data = pd.read_csv(file_path, sep="\t")

results = []

for index, row in data.iterrows():
    chrom = row[1]
    start = row[2]
    end = row[3]

    print(f"\nОбработка региона: {chrom}:{start}-{end}")

    nearest_transcript_info = find_nearest_transcript(chrom, start, end)
    if nearest_transcript_info is not None:
        nearest_transcript, strand, biotype = nearest_transcript_info
        print("Ближайший транскрипт:")
        print(f"ID: {nearest_transcript.transcript_id}")
        print(f"Ген: {nearest_transcript.gene_name}")
        print(f"Биотип: {biotype}")
        print(f"Цепь: {strand}")

        results.append([chrom, start, end, nearest_transcript.transcript_id, 
                       nearest_transcript.gene_name, strand, biotype])
    else:
        print("Транскрипт не найден.")
        # регион без транскрипта
        results.append([chrom, start, end, "N/A", "N/A", "N/A", "N/A"])

results_df = pd.DataFrame(results, columns=["chrom", "start", "end", "transcript_id", 
                                           "gene_name", "strand", "biotype"])
output_file = "icrs_nearest_tr.csv"
results_df.to_csv(output_file, index=False)

print(f"Результаты сохранены в файл: {output_file}")
