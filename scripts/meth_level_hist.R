library(data.table)
library(ggplot2)
library(R.utils)

file_path <- "merged_hg38.bed"

data <- fread(file_path, header = FALSE)
data_m <- subset(data, V4 == 'm')
data_m <- subset(data_m, V1 == c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX','chrM'))

methylation_levels <- as.numeric(data_m$V11) / 100
methylation_levels <- methylation_levels[!is.na(methylation_levels) & methylation_levels >= 0 & methylation_levels <= 1]

png("autosome", width = 800, height = 600)

p <- hist(methylation_levels, 
        breaks = seq(0, 1, by = 0.02), freq = TRUE,
        col = "gray", 
        border = "black", 
        main = "DNA methylation in autosomes", 
        xlab = "DNA Methylation level")

dev.off()
print(p)

