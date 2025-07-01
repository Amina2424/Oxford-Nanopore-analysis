library(data.table)
library(ggplot2)
library(R.utils)


file_path <- "haplotype_ungrouped.bed"

data <- fread(file_path, header = FALSE)

methylation_levels <- as.numeric(data$V11) / 100

methylation_levels <- methylation_levels[!is.na(methylation_levels) & methylation_levels >= 0 & methylation_levels <= 1]

ggplot(data.frame(methylation_levels), aes(x = methylation_levels)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "Methylation Level Distribution", 
       x = "Methylation Level", 
       y = "Count") +
  theme_minimal()
