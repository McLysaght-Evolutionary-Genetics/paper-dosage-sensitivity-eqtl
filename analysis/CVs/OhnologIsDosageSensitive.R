library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstatix)
library(stringr)

#geneList
singletons <- read_tsv("../datasets/geneLists/singletons.txt")
ohnologs <- read_tsv("../datasets/geneLists/Singh and Isambert/hsapiens.Pairs.Relaxed.2R.Ens75.1-Y.txt")
SSds <- read_tsv("../datasets/geneLists/SSDsRelaxedOhnos.txt")

selected_gene <- rbind(
  singletons |> mutate(type = "singleton"),
  ohnologs |> mutate(type = "ohnolog"), 
  SSds |> mutate(type = "SSD")
)

#tissue list

tissues <- list.files("../datasets/GTEx_Analysis_v7_eQTL_expression_matrices/", pattern = "bed\\.gz$") |> 
  str_remove(".v7.normalized_expression.bed.gz")

##TPM from GTEX

samples <- read_tsv("GTEx_v7_Annotations_SampleAttributesDS.txt")
library(data.table)
TPM <- fread("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz", skip = 2)

filtered_TPM <- TPM |> 
  mutate(Name = str_remove(Name, "\\..*$")) |> 
  inner_join(selected_gene, by = join_by(Name == `Ensembl Gene ID`))

filtered_TPM <- filtered_TPM |> 
  pivot_longer(cols = -c(Name, Description, type), names_to = "SAMPID") |> 
  inner_join(select(samples, SAMPID, SMTSD))

filtered_TPM <- filtered_TPM %>%
  mutate(log_value = log2(value + 1))


#calculate CV
cv_by_gene_tissue <- filtered_TPM %>%
  group_by(Name, SMTSD) %>%
  summarise(
    mean_expr = mean(log_value, na.rm = TRUE),
    sd_expr   = sd(log_value, na.rm = TRUE),
    cv_expr   = ifelse(mean_expr > 0, sd_expr / mean_expr, NA_real_),
    .groups = "drop"
  )


# tissues

sample_counts <- filtered_TPM |> 
  select(SAMPID, SMTSD) |> 
  distinct() |> 
  group_by(SMTSD) |> 
  summarise(n = n())

sample_counts <- sample_counts |> 
  mutate(tissue_name = str_replace_all(SMTSD, "[^A-Za-z0-9]+", "_") |> 
           str_remove("_+$") |> 
           str_replace("cervical_c_1", "cervical_c-1") |> 
           str_replace("EBV_transformed", "EBV-transformed")) |> 
  filter(tissue_name %in% tissues)


# plot 1, ohnologs

cv_clean <- cv_by_gene_tissue |> 
  filter(is.finite(cv_expr)) |> 
  inner_join(selected_gene, by = join_by(Name == `Ensembl Gene ID`))
  
cv_labeled <- cv_clean |> 
  inner_join(sample_counts, by = "SMTSD") |> 
  mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")"))


wilcox_results <- cv_clean |> 
  group_by(SMTSD) |> 
  wilcox_test(cv_expr ~ type, 
              comparisons = list(
                c("ohnolog", "singleton"),
                c("ohnolog", "SSD"),
                c("SSD", "singleton")
              )) |> 
  adjust_pvalue(method = "BH") |> 
  mutate(p.signif = p_format(p, digits = 2))

wilcox_results <- wilcox_results |> 
  mutate(y.position = rep(c(1.1, 1.25, 1.4), times = length(unique(SMTSD))))

cv_labeled |> 
  ggplot(aes(x = type, y = cv_expr)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~SMTSD_label, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    data = wilcox_results |> 
      inner_join(sample_counts, by = "SMTSD") |> 
      mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")")),
    label = "p.adj.signif",
    y.position = "y.position"
  ) +
  coord_cartesian(ylim = c(0, 1.5)) +
  xlab("") + 
  ylab("CV（Coefficient of Variation）") + 
  theme_cowplot()

ggsave("cv_ohnolog.pdf", width = 18, height = 27)

#plot2 CNV
CNV <- read_tsv("../outputFiles/genesWitheQTLTissueCountMetasoftAndCNVExACStatus.csv") |> 
  select(2, CNV)

cv_CNV_clean <- cv_by_gene_tissue |> 
  filter(is.finite(cv_expr)) |> 
  inner_join(CNV, by = join_by(Name == `Ensembl Gene ID`))

cv_CNV_labeled <- cv_CNV_clean |> 
  inner_join(sample_counts, by = "SMTSD") |> 
  mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")"))


wilcox_CNV_results <- cv_CNV_clean |> 
  group_by(SMTSD) |> 
  wilcox_test(cv_expr ~ CNV, 
              comparisons = list(
                c("N", "Y")
              )) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance("p.adj") |> 
  mutate(p.signif = p_format(p, digits = 2))

wilcox_CNV_results <- wilcox_CNV_results |> 
  mutate(y.position = rep(c(0.7), times = length(unique(SMTSD))))

cv_CNV_labeled |> 
  ggplot(aes(x = CNV, y = cv_expr)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~SMTSD_label, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    data = wilcox_CNV_results |> 
      inner_join(sample_counts, by = "SMTSD") |> 
      mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")")),
    label = "p.adj.signif",
    y.position = "y.position"
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  xlab("") + 
  ylab("CV（Coefficient of Variation）") +
  scale_x_discrete(labels = c("N" = "CNV-free", "Y" = "CNV-affected")) +
  theme_cowplot()

ggsave("cv_CNV.pdf", width = 18, height = 27)

#plot3 Haploinsufficient

haplo <- read_tsv("../outputFiles/genesWitheQTLTissueCountMetasoftAndHaploStatus.csv") |> 
  select(2, haplo)

cv_haplo_clean <- cv_by_gene_tissue |> 
  filter(is.finite(cv_expr)) |> 
  inner_join(haplo, by = join_by(Name == `Ensembl Gene ID`))

cv_haplo_labeled <- cv_haplo_clean |> 
  inner_join(sample_counts, by = "SMTSD") |> 
  mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")"))


wilcox_haplo_results <- cv_haplo_clean |> 
  group_by(SMTSD) |> 
  wilcox_test(cv_expr ~ haplo, 
              comparisons = list(
                c("N", "Y")
              )) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance("p.adj") |> 
  mutate(p.signif = p_format(p, digits = 2))

wilcox_haplo_results <- wilcox_haplo_results |> 
  mutate(y.position = rep(c(0.7), times = length(unique(SMTSD))))

cv_haplo_labeled |> 
  ggplot(aes(x = haplo, y = cv_expr)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~SMTSD_label, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    data = wilcox_haplo_results |> 
      inner_join(sample_counts, by = "SMTSD") |> 
      mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")")),
    label = "p.adj.signif",
    y.position = "y.position"
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  xlab("") + 
  ylab("CV（Coefficient of Variation）") +
  scale_x_discrete(labels = c("N" = "others", "Y" = "Haplo-\ninsufficient")) +
  theme_cowplot()

ggsave("cv_haplo.pdf", width = 18, height = 27)

#plot4 CCN genes

CCN <- read_tsv("../outputFiles/genesWitheQTLTissueCountMetasoftAndCCNStatus.csv") |> 
  select(2, CCN)

cv_CCN_clean <- cv_by_gene_tissue |> 
  filter(is.finite(cv_expr)) |> 
  inner_join(CCN, by = join_by(Name == `Ensembl Gene ID`))

cv_CCN_labeled <- cv_CCN_clean |> 
  inner_join(sample_counts, by = "SMTSD") |> 
  mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")"))


wilcox_CCN_results <- cv_CCN_clean |> 
  group_by(SMTSD) |> 
  wilcox_test(cv_expr ~ CCN, 
              comparisons = list(
                c("N", "Y")
              )) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance("p.adj") |> 
  mutate(p.signif = p_format(p, digits = 2))

wilcox_CCN_results <- wilcox_CCN_results |> 
  mutate(y.position = rep(c(0.7), times = length(unique(SMTSD))))

cv_CCN_labeled |> 
  ggplot(aes(x = CCN, y = cv_expr)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~SMTSD_label, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    data = wilcox_CCN_results |> 
      inner_join(sample_counts, by = "SMTSD") |> 
      mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")")),
    label = "p.adj.signif",
    y.position = "y.position"
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  xlab("") + 
  ylab("CV（Coefficient of Variation）") +
  scale_x_discrete(labels = c("N" = "not conserved", "Y" = "CCN")) +
  theme_cowplot()

ggsave("cv_CCN.pdf", width = 18, height = 27)

#plot5 CNVRs genes

CNVRs <- read_tsv("../outputFiles/genesWitheQTLTissueCountMetasoftAndCNVZarreiStatus.csv") |> 
  select(2, CNV) |> 
  rename(CNVRs = CNV)

cv_CNVRs_clean <- cv_by_gene_tissue |> 
  filter(is.finite(cv_expr)) |> 
  inner_join(CNVRs, by = join_by(Name == `Ensembl Gene ID`))

cv_CNVRs_labeled <- cv_CNVRs_clean |> 
  inner_join(sample_counts, by = "SMTSD") |> 
  mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")"))


wilcox_CNVRs_results <- cv_CNVRs_clean |> 
  group_by(SMTSD) |> 
  wilcox_test(cv_expr ~ CNVRs, 
              comparisons = list(
                c("N", "Y")
              )) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance("p.adj") |> 
  mutate(p.signif = p_format(p, digits = 2))

wilcox_CNVRs_results <- wilcox_CNVRs_results |> 
  mutate(y.position = rep(c(0.7), times = length(unique(SMTSD))))

cv_CNVRs_labeled |> 
  ggplot(aes(x = CNVRs, y = cv_expr)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~SMTSD_label, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    data = wilcox_CNVRs_results |> 
      inner_join(sample_counts, by = "SMTSD") |> 
      mutate(SMTSD_label = paste0(SMTSD, " \n(n = ", n, ")")),
    label = "p.adj.signif",
    y.position = "y.position"
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  xlab("") + 
  ylab("CV（Coefficient of Variation）") +
  scale_x_discrete(labels = c("N" = "outside CNVRs", "Y" = "CNVR genes")) +
  theme_cowplot()

ggsave("cv_CNVRs.pdf", width = 18, height = 27)


# save  -------------------------------------------------------------------
library(openxlsx)

list(wilcox_results = wilcox_results, 
     wilcox_haplo_results = wilcox_haplo_results, 
     wilcox_CNVRs_results = wilcox_CNVRs_results, 
     wilcox_CNV_results = wilcox_CNV_results, 
     wilcox_CCN_results = wilcox_CCN_results) |> 
  write.xlsx(file = "wilcox_results.xlsx", rowNames = FALSE)
