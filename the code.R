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

coldata <- read.csv("coldata.csv")



rownames(counts) <- counts[, 1]
counts <- counts[, -1]

rownames(coldata) <- coldata[, 1]
coldata <- coldata[, -1]


### comparing WT vs MCT1 KO ANG/PE

coldata_drug <- coldata[-(1:17), ]
counts_drug <- counts[ ,-(1:17)]

all(colnames(counts) == rownames(coldata))

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

# Comparing WT Saline vs. WT ANG/PE

to_delete <- c(1:7, 14:17, 21:25)

coldata_WT <- coldata[-to_delete,]
counts_WT <- counts[,-to_delete]

all(colnames(counts) == rownames(coldata))

dds_WT <- DESeqDataSetFromMatrix(countData = counts_WT,
                              colData = coldata_WT,
                              design = ~ Infusion)

smallestGroupSize <- 3
keep <- rowSums(counts(dds_WT) >= 10) >= smallestGroupSize
dds_WT <- dds_WT[keep,]

dds_KO$Infusion <- factor(dds_KO$Infusion, levels = c("Saline","ANGPE"))

dds_WT <- DESeq(dds_WT)

res_WT <- results(dds_WT)

res_df_WT <- as.data.frame(res_WT)

View(res_df_WT)

write.csv(res_df_WT, "WT_SalinevsANGPE_unfiltered_flipped.csv")

res_df$log10_padj <- -log10(res_df$padj)

res_df |>
  ggplot(aes(x = log2FoldChange, y = log10_padj)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 padj") +
  theme_minimal() +
  geom_point()

# Comparing WT Saline to KO Saline

to_delete <- c(1:7, 18:25)

coldata_Saline <- coldata[-to_delete,]
counts_Saline <- counts[,-to_delete]

all(colnames(counts_Saline) == rownames(coldata_Saline))

dds_Saline <- DESeqDataSetFromMatrix(countData = counts_Saline,
                                 colData = coldata_Saline,
                                 design = ~ Genotype)

smallestGroupSize <- 3
keep <- rowSums(counts(dds_Saline) >= 10) >= smallestGroupSize
dds_Saline <- dds_Saline[keep,]

dds_Saline$Genotype <- factor(dds_Saline$Genotype, levels = c("WT","MCT1 KO"))

dds_Saline <- DESeq(dds_Saline)

res_Saline <- results(dds_Saline)

res_df_Saline <- as.data.frame(res_Saline)

View(res_df_Saline)

write.csv(res_df_Saline, "WTSalinevsKOSaline_filtered.csv")

# Comparing KO Saline to KO Drug

to_delete <- c(1:13, 18:20)

coldata_KO <- coldata[-to_delete,]
counts_KO <- counts[,-to_delete]

all(colnames(counts_KO) == rownames(coldata_KO))

dds_KO <- DESeqDataSetFromMatrix(countData = counts_KO,
                                     colData = coldata_KO,
                                     design = ~ Infusion)

smallestGroupSize <- 3
keep <- rowSums(counts(dds_KO) >= 10) >= smallestGroupSize
dds_KO <- dds_KO[keep,]

dds_KO$Infusion <- factor(dds_KO$Infusion, levels = c("Saline","ANGPE"))

dds_KO <- DESeq(dds_KO)

res_KO <- results(dds_KO)

res_df_KO <- as.data.frame(res_KO)

View(res_df_KO)

write.csv(res_df_KO, "KOSalinevsKODrug_filtered.csv")






