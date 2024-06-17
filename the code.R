if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

library(tidyverse)
library(DESeq2)
library(purrr)

# Import Inefficiently ----------------------------------------------------

one <- read_tsv("1_WTNOInfusion.tsv")

two <- read_tsv("2_WTNoInfusion.tsv")

three <- read_tsv("3_WTNoInfusion.tsv")

four <- read_tsv("4_WTNoInfusion.tsv")

five <- read_tsv("5_MCT1KONoInfusion.tsv")

six <- read_tsv("6_MCT1KONoInfusion.tsv")

seven <- read_tsv("7_MCT1KONoInfusion.tsv")

eight <- read_tsv("8_WTSalineInfusion.tsv")

nine <- read_tsv("9_WTSalineInfusion.tsv")

ten <- read_tsv("10_WTSalineInfusion.tsv")

eleven <- read_tsv("11_WTSalineInfusion.tsv")

twelve <- read_tsv("12_WTSalineInfusion.tsv")

thirteen <- read_tsv("13_WTSalineInfusion.tsv")

fourteen <- read_tsv("14_MCT1KOSalineInfusion.tsv")

fifteen <- read_tsv("15_MCT1KOSalineInfusion.tsv")

sixteen <- read_tsv("16_MCT1KOSalineInfusion.tsv")

seventeen <- read_tsv("17_MCT1KOSalineInfusion.tsv")

eighteen <- read_tsv("18_WTANGPEInfusion.tsv")

nineteen <- read_tsv("19_WTANGPEInfusion.tsv")

twenty <- read_tsv("20_WTANGPEInfusion.tsv")

twentyone <- read_tsv("21_MCT1KOANGPEInfusion.tsv")

twentytwo <- read_tsv("22_MCT1KOANGPEInfusion.tsv")

twentythree <- read_tsv("23_MCT1KOANGPEInfusion.tsv")

twentyfour <- read_tsv("24_MCT1KOANGPEInfusion.tsv")

twentyfive <- read_tsv("25_MCT1KOANGPEInfusion.tsv")


one <- one |>
  rename_at(vars(2), ~ "1_WTNoInfusion")

two <- two |>
  rename_at(vars(2), ~ "2_WTNoInfusion")

three <- three |>
  rename_at(vars(2), ~ "3_WTNoInfusion")

four <- four |>
  rename_at(vars(2), ~ "4_WTNoInfusion")

five <- five |>
  rename_at(vars(2), ~ "5_MCT1KONoInfusion")

six <- six |>
  rename_at(vars(2), ~ "6_MCT1KONoInfusion")

seven <- seven |>
  rename_at(vars(2), ~ "7_MCT1KONoInfusion")

eight <- eight |>
  rename_at(vars(2), ~ "8_WTSalineInfusion")

nine <- nine |>
  rename_at(vars(2), ~ "9_WTSalineInfusion")

ten <- ten |>
  rename_at(vars(2), ~ "10_WTSalineInfusion")

eleven <- eleven |>
  rename_at(vars(2), ~ "11_WTSalineInfusion")

twelve <- twelve |>
  rename_at(vars(2), ~ "12_WTSalineInfusion")

thirteen <- thirteen |>
  rename_at(vars(2), ~ "13_WTSalineInfusion")

fourteen <- fourteen |>
  rename_at(vars(2), ~ "14_MCT1KOSalineInfusion")

fifteen <- fifteen |>
  rename_at(vars(2), ~ "15_MCT1KOSalineInfusion")

sixteen <- sixteen |>
  rename_at(vars(2), ~ "16_MCT1KOSalineInfusion")

seventeen <- seventeen |>
  rename_at(vars(2), ~ "17_MCT1KOSalineInfusion")

eighteen <- eighteen |>
  rename_at(vars(2), ~ "18_WTANGPEInfusion")

nineteen <- nineteen |>
  rename_at(vars(2), ~ "19_WTANGPEInfusion")

twenty <- twenty |>
  rename_at(vars(2), ~ "20_WTANGPEInfusion")

twentyone <- twentyone |>
  rename_at(vars(2), ~ "21_MCT1KOANGPEInfusion")

twentytwo <- twentytwo |>
  rename_at(vars(2), ~ "22_MCT1KOANGPEInfusion")

twentythree <- twentythree |>
  rename_at(vars(2), ~ "23_MCT1KOANGPEInfusion")

twentyfour <- twentyfour |>
  rename_at(vars(2), ~ "24_MCT1KOANGPEInfusion")

twentyfive <- twentyfive |>
  rename_at(vars(2), ~ "25_MCT1KOANGPEInfusion")

df_list <- list(one, two, three, four, five, six, seven, eight, nine, ten,
                eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen, nineteen, twenty,
                twentyone, twentytwo, twentythree, twentyfour, twentyfive)


combined_df <- df_list[[1]]
for (i in 2:length(df_list)) {
  combined_df <- left_join(combined_df, df_list[[i]], by = "Geneid")
}

counts <- as.data.frame(combined_df)

rownames(counts) <- counts[, 1]
counts <- counts[, -1]


coldata <- read.csv("coldata.csv")

rownames(coldata) <- coldata[, 1]
coldata <- coldata[, -1]

coldata_drug <- coldata[-(1:17), ]
counts_drug <- counts[ ,-(1:17)]

all(colnames(counts) == rownames(coldata))

### comparing WT vs MCT1 KO ANG/PE

dds <- DESeqDataSetFromMatrix(countData = counts_drug,
                              colData = coldata_drug,
                              design = ~ Genotype)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$Genotype <- factor(dds$Genotype, levels = c("WT","MCT1 KO"))


dds <- DESeq(dds)

res <- results(dds)

res_df <- as.data.frame(res)

View(res_df)

write.csv(res_df, "ANGPE_filtered.csv")

res_df$log10_padj <- -log10(res_df$padj)

res_df |>
  ggplot(aes(x = log2FoldChange, y = log10_padj)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 padj") +
  theme_minimal() +
  geom_point()

### comparing WT vs MCT1 KO ANG/PE - No filter - probably use this

dds <- DESeqDataSetFromMatrix(countData = counts_drug,
                              colData = coldata_drug,
                              design = ~ Genotype)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$Genotype <- factor(dds$Genotype, levels = c("WT","MCT1 KO"))


dds <- DESeq(dds)

res <- results(dds, independentFiltering=FALSE)

res_df <- as.data.frame(res)

View(res_df)

write.csv(res_df, "ANGPE_unfiltered.csv")

res_df$log10_padj <- -log10(res_df$padj)

res_df |>
  ggplot(aes(x = log2FoldChange, y = log10_padj)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 padj") +
  theme_minimal() +
  geom_point()




