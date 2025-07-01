table <- read.delim('DMRs_330cov.csv')
table <- table[-14,]
boxplot(table$percent_m, xlab = "Percentation of DMRs methylation", ylab = "Percent (%)",
        col.axis = "darkgreen", col.lab = "darkgreen", border = 'darkgreen', col = 'lightgreen')

sd <- sd(table$percent_m) # [1] 5.157023
mean <- mean(table$percent_m) #[1] 55.80415

data <- table[order(table$percent_m, decreasing = TRUE), ]



bar_positions <- barplot(height = data$percent_m,
                         names.arg = data$name,
                         las = 2,  # Повернем подписи на 90 градусов
                         col = "lightgreen",
                         main = "Methylation of Differential Methylated Regions (gDMRs)",
                         ylab = " Percent (%)",
                         cex.names = 0.45)  

text(x = bar_positions, 
     y = data$percent_m, 
     labels = paste0(round(data$percent_m), "%"), 
     pos = 1, 
     cex = 0.5, 
     offset = 3,
     col = "darkgreen")


table$name <- as.factor(table$name)

ggplot(table, aes(x = name, y = percent_m, fill = name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.5, fill = "white") +
  labs(title = "DMRs' Methylation Percent ",
       x = "Region",
       y = "Methylation Percent") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
        legend.position = "none")


ggplot(table, aes(x = name, y = percent_m) +
  geom_violin() +  
  labs(title = "Метилирование для регионов",
       x = "Регион",
       y = "Процент метилирования") +
  theme_minimal()  +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none" 
  )

ggplot(df, aes(x = factor(name), y = percent_m)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") + 
  labs(title = "Violin Plot for name and percent_m",
        x = "Name",
       y = "Percent Methylation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  