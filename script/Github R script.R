####################################
# Indoor Airborne Microbiome and Lung Health Study
####################################

#---------------------------------
# Set your working directory
#---------------------------------
# Load required libraries
#---------------------------------
library(phyloseq)
library(ggplot2)
library(readxl)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(microbiome)
library(data.table)
library(table1)
library(htmltools)
library(ampvis2)
library(broom)
library(compositions)
library(writexl)
library(vegan)

#---------------------------------
# Construct Phyloseq object
#---------------------------------
otu_table_data <- fread("otu_table.csv", data.table = FALSE)
rownames(otu_table_data) <- otu_table_data[,1]
otu_table_data <- otu_table_data[,-1]

tax_table_data <- read.csv("tax_table.csv", row.names = 1)
sample_data_data <- read.csv("sample_metadata.csv", sep = ',', header = TRUE, row.names = 1)

otu_tab <- otu_table(as.matrix(otu_table_data), taxa_are_rows = FALSE)
tax_tab <- tax_table(as.matrix(tax_table_data))
sample_dat <- sample_data(sample_data_data)

ps <- phyloseq(otu_tab, tax_tab, sample_dat)

#---------------------------------
# Table 1 / Participant Characteristics
#---------------------------------
lung_data <- read_excel("sample_metadata.xlsx")

lung_data$gender <- factor(lung_data$gender, c("male", "female", ""))
lung_data$atopy <- factor(lung_data$atopy, c("atopic", "non_atopic", ""))
lung_data$smoker <- factor(lung_data$smoker, c("never_smoker", "ex_smoker", "current_smoker"))

label(lung_data$gender) <- "Gender"
label(lung_data$smoker) <- "Smoking"
label(lung_data$feno) <- "FeNO (ppb)"
label(lung_data$weight) <- "Weight"
label(lung_data$height) <- "Height"

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2),
       c("", "Mean (SD)" = sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

charcter.table.1 <- table1(
  ~ age + weight + height + smoker + Z_score.FEV1 + Z_score.FVC + Z_score.FEV1FVC + feno + atopy | gender,
  data = lung_data,
  render.continuous = my.render.cont,
  render.categorical = my.render.cat,
  overall = "Total",
  topclass = "Rtable1-grid"
)

write.csv(charcter.table.1, "Table 1.csv", row.names = FALSE)

#---------------------------------
# Figure 1 / Taxonomic Heatmaps
#---------------------------------
otutable <- data.frame(
  OTU = colnames(phyloseq::otu_table(ps)@.Data),
  t(phyloseq::otu_table(ps)@.Data),
  phyloseq::tax_table(ps)@.Data,
  check.names = FALSE
)
row.names(otutable) <- paste("ASV", 1:dim(otutable)[1], sep = "")
otutable$OTU <- row.names(otutable)

metadata <- read_excel("sample_metadata.xlsx")
ampvis <- amp_load(otutable = otutable, metadata = metadata)

p1 <- amp_heatmap(ampvis, group_by = "ethnicity", tax_aggregate = "Phylum", tax_show = 10,
                  showRemainingTaxa = TRUE, plot_colorscale = "sqrt", plot_values = TRUE,
                  min_abundance = 0.1, measure = "mean") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=12)) +
  theme(legend.position = "none")

p2 <- amp_heatmap(ampvis, group_by = "ethnicity", tax_aggregate = "Class", tax_show = 15,
                  showRemainingTaxa = TRUE, plot_colorscale = "sqrt", plot_values = TRUE,
                  min_abundance = 0.1, measure = "mean") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=12)) +
  theme(legend.position = "none")

p3 <- amp_heatmap(ampvis, group_by = "ethnicity", tax_aggregate = "Genus", tax_show = 30,
                  showRemainingTaxa = TRUE, plot_colorscale = "sqrt", plot_values = TRUE,
                  min_abundance = 0.1, measure = "mean") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=11)) +
  theme(legend.position = "right")

heatmap <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(2,2,2), labels = c("A", "B", "C"))

ggsave("Figure 1.jpeg", plot = heatmap, height = 15, width = 15, units = 'in', dpi = 300)

#---------------------------------
# Table 2 / Association Between Microbial Diversity, Load, and Outcomes
#---------------------------------
exposures <- c("diversity_shannon", "rcopies_m2", "eu_m2")
outcomes <- c("Z_score.FEV1", "Z_score.FVC", "Z_score.FEV1FVC", "feno")

for (exp in exposures) {
  iqr_value <- IQR(lung_data[[exp]], na.rm = TRUE)
  lung_data[[paste0(exp, "_iqr")]] <- lung_data[[exp]] / iqr_value
}

lung_data_male <- subset(lung_data, gender == "male")
lung_data_female <- subset(lung_data, gender == "female")

models <- list()

for (exp in exposures) {
  exp_var <- paste0(exp, "_iqr")
  
  for (outcome in outcomes) {
    
    if (outcome == "feno") {
      formula_male <- as.formula(paste(outcome, "~", exp_var, "+ atopy + sample_center + age + weight + smoker"))
      formula_female <- as.formula(paste(outcome, "~", exp_var, "+ atopy + sample_center + age + weight + smoker"))
    } else {
      formula_male <- as.formula(paste(outcome, "~", exp_var, "+ atopy + sample_center + weight + smoker"))
      formula_female <- as.formula(paste(outcome, "~", exp_var, "+ atopy + sample_center + weight + smoker"))
    }
    
    models[[paste0(outcome, "_male_", exp)]] <- lm(formula_male, data = lung_data_male)
    models[[paste0(outcome, "_female_", exp)]] <- lm(formula_female, data = lung_data_female)
  }
}

extract_results <- function(model) {
  coefs <- summary(model)$coefficients
  beta <- coefs[2, "Estimate"]
  lower <- confint(model)[2, 1]
  upper <- confint(model)[2, 2]
  pval <- coefs[2, "Pr(>|t|)"]
  
  data.frame(
    Beta = round(beta, 2),
    CI_lower = round(lower, 2),
    CI_upper = round(upper, 2),
    P_value = round(pval, 3)
  )
}

results_list <- list()

for (exp in exposures) {
  for (outcome in outcomes) {
    male_model <- models[[paste0(outcome, "_male_", exp)]]
    female_model <- models[[paste0(outcome, "_female_", exp)]]
    
    male_result <- extract_results(male_model)
    female_result <- extract_results(female_model)
    
    results_list[[paste0(exp, "_", outcome)]] <- cbind(
      Exposure = exp,
      Outcome = outcome,
      Sex = c("Men", "Women"),
      rbind(male_result, female_result)
    )
  }
}

final_results <- do.call(rbind, results_list)

write.csv(final_results, "Table 2.csv", row.names = FALSE)

#---------------------------------
# Figure 2 / Association Between Richness Within Bacterial Classes and Outcomes
#---------------------------------

# Define classes and outcomes
classes <- c("Alphaproteobacteria", "Actinobacteria", "Bacilli", 
             "Bacteroidia", "Clostridia", "Gammaproteobacteria")

outcomes <- c("Z_score.FEV1", "Z_score.FVC", "Z_score.FEV1FVC", "feno")

# Prepare empty list
results_list <- list()

# Loop over classes
for (cls in classes) {
  ps_sub <- subset_taxa(ps, Class == cls)
  richness_df <- alpha(ps_sub, index = "observed")
  richness_df$sam_name <- rownames(richness_df)
  
  meta_df <- meta(ps)
  meta_df$sam_name <- rownames(meta_df)
  
  df <- merge(richness_df, meta_df, by = "sam_name")
  df$Class <- cls
  
  for (sex in c("male", "female")) {
    df_sex <- subset(df, gender == sex)
    
    for (outcome in outcomes) {
      if (outcome == "feno") {
        formula <- as.formula("feno ~ observed + atopy + sample_center + age + weight + smoker")
      } else {
        formula <- as.formula(paste0(outcome, " ~ observed + atopy + sample_center + weight + smoker"))
      }
      
      model <- lm(formula, data = df_sex)
      tidy_out <- tidy(model)[2, ]
      
      tidy_out$Outcome <- ifelse(outcome == "feno", "FeNO", gsub("Z_score.", "", outcome))
      tidy_out$Sex <- ifelse(sex == "male", "(men)", "(women)")
      tidy_out$Outcome <- paste0(tidy_out$Outcome, " ", tidy_out$Sex)
      tidy_out$Class <- cls
      
      results_list[[length(results_list) + 1]] <- tidy_out
    }
  }
}

# Combine all results
figure2_df <- bind_rows(results_list)

# Create Direction column
figure2_df$Direction <- ifelse(figure2_df$p.value < 0.05 & figure2_df$estimate > 0, "Positive",
                               ifelse(figure2_df$p.value < 0.05 & figure2_df$estimate < 0, "Negative", "None"))

# Set factor levels
figure2_df$Class <- factor(figure2_df$Class, levels = classes)
figure2_df$Outcome <- factor(figure2_df$Outcome, levels = c(
  "FEV1 (men)", "FEV1 (women)", 
  "FVC (men)", "FVC (women)", 
  "FEV1FVC (men)", "FEV1FVC (women)", 
  "FeNO (men)", "FeNO (women)"
))

# Plot Figure 2
Figure.2 <- ggplot(figure2_df, aes(x = Outcome, y = Class, fill = Direction)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_manual(
    values = c("Positive" = "#0072B2", "Negative" = "#E69F00", "None" = "white"),
    name = "Association"
  ) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_x_discrete(
    labels = c(
      "FEV1 (men)" = expression(FEV[1]*" (men)"),
      "FEV1 (women)" = expression(FEV[1]*" (women)"),
      "FVC (men)" = "FVC (men)",
      "FVC (women)" = "FVC (women)",
      "FEV1FVC (men)" = expression(FEV[1]*"/FVC (men)"),
      "FEV1FVC (women)" = expression(FEV[1]*"/FVC (women)"),
      "FeNO (men)" = "FeNO (men)",
      "FeNO (women)" = "FeNO (women)"
    )
  )

ggsave(filename = "Figure 2.jpeg", plot = Figure.2, height = 7, width = 9, units = 'in', dpi = 300)

write.csv(figure2_df, "bacterial_richness_respiratory_outcomes.csv", row.names = FALSE)

#---------------------------------
# Table 3 / Mantel Test Between Beta-Diversity and Outcomes
#---------------------------------

mantel_results <- list()

for (outcome in outcomes) {
  ps_sub <- subset_samples(ps, !is.na(get_variable(ps, outcome)))
  
  # Male
  ps_male <- subset_samples(ps_sub, sample_data(ps_sub)$gender == "male")
  otu_male <- otu_table(ps_male)
  high_otu_idx <- colSums(otu_male)/sum(otu_male) >= 0.00001
  ps_male <- prune_taxa(taxa_names(ps_male)[high_otu_idx], ps_male)
  
  a_male <- otu_table(ps_male)
  b_male <- a_male[, colSums(a_male) > 0] + 0.001
  beta_male <- coda.base::dist(b_male, method = "aitchison")
  
  outcome_male <- sample_data(ps_male)[[outcome]]
  dist_male <- dist(outcome_male, method = "euclidean")
  
  mantel_male <- mantel(beta_male, dist_male, method = "spearman", permutations = 9999, na.rm = TRUE)
  
  # Female
  ps_female <- subset_samples(ps_sub, sample_data(ps_sub)$gender == "female")
  otu_female <- otu_table(ps_female)
  high_otu_idx <- colSums(otu_female)/sum(otu_female) >= 0.00001
  ps_female <- prune_taxa(taxa_names(ps_female)[high_otu_idx], ps_female)
  
  a_female <- otu_table(ps_female)
  b_female <- a_female[, colSums(a_female) > 0] + 0.001
  beta_female <- coda.base::dist(b_female, method = "aitchison")
  
  outcome_female <- sample_data(ps_female)[[outcome]]
  dist_female <- dist(outcome_female, method = "euclidean")
  
  mantel_female <- mantel(beta_female, dist_female, method = "spearman", permutations = 9999, na.rm = TRUE)
  
  # Collect
  mantel_results[[outcome]] <- data.frame(
    Outcome = outcome,
    R_male = round(mantel_male$statistic, 3),
    P_male = signif(mantel_male$signif, 2),
    R_female = round(mantel_female$statistic, 3),
    P_female = signif(mantel_female$signif, 2)
  )
}

mantel_summary <- bind_rows(mantel_results)

write.csv(mantel_summary, "Table 3.csv", row.names = FALSE)

#---------------------------------
# Figure 3 / Genera Associations with Outcomes
#---------------------------------

metadata_all <- read_excel("sample_metadata.xlsx")
relative_abundance <- read_excel("rel_abund.genera.dust.xlsx")

# Taxonomy info
tax_table_df <- as.data.frame(tax_table(ps))
tax_table_df$Genus <- rownames(tax_table_df)
taxonomy_info <- tax_table_df %>% select(Genus, Phylum, Class)

# Split by sex
data.male <- filter(relative_abundance, gender == "male")
data.female <- filter(relative_abundance, gender == "female")

# Define CLR + IQR function
run_clr_lm <- function(data, metadata, outcome_var, adjust_vars) {
  df <- inner_join(metadata, data, by = "sample_id")
  
  features_only <- df %>% select(where(is.numeric)) %>% select(-all_of(outcome_var))
  
  clr_transformed <- clr(features_only + 1e-6) %>% as.data.frame()
  scaled_clr <- scale(clr_transformed, center = FALSE, scale = apply(clr_transformed, 2, IQR)) %>% as.data.frame()
  
  full_data <- cbind(
    df %>% select(all_of(outcome_var), all_of(adjust_vars)),
    scaled_clr
  )
  
  results <- lapply(names(scaled_clr), function(feature) {
    formula <- as.formula(paste(outcome_var, "~", feature, "+", paste(adjust_vars, collapse = "+")))
    model <- lm(formula, data = full_data)
    tidy(model)[2, ]
  }) %>% bind_rows()
  
  results$Genus <- names(scaled_clr)
  return(results)
}

# Run models
all_results <- list()

for (analysis in analyses) {
  outcome_var <- analysis$outcome
  adj <- if (!is.null(analysis$adjust)) analysis$adjust else adjust_vars
  
  metadata.selected <- metadata_all %>%
    select(sample_id, all_of(outcome_var), all_of(adj))
  
  metadata.male <- metadata.selected %>%
    filter(sample_id %in% data.male$sample_id) %>%
    mutate(across(all_of(adj), as.factor))
  
  metadata.female <- metadata.selected %>%
    filter(sample_id %in% data.female$sample_id) %>%
    mutate(across(all_of(adj), as.factor))
  
  res_male <- run_clr_lm(data.male, metadata.male, outcome_var, adj) %>%
    filter(p.value < 0.05) %>%
    mutate(Outcome = analysis$label, Sex = "men")
  
  res_female <- run_clr_lm(data.female, metadata.female, outcome_var, adj) %>%
    filter(p.value < 0.05) %>%
    mutate(Outcome = analysis$label, Sex = "women")
  
  all_results[[length(all_results)+1]] <- bind_rows(res_male, res_female)
}

# Combine
final_df <- bind_rows(all_results)

# Add Direction column
final_df <- final_df %>%
  mutate(Direction = case_when(
    estimate > 0 ~ "Positive",
    estimate < 0 ~ "Negative",
    TRUE ~ NA_character_
  ))

final_df <- final_df %>%
  mutate(Outcome_Sex = paste0(Outcome, " (", Sex, ")"))

# Prepare for heatmap
all_combinations <- expand.grid(
  Genus = unique(final_df$Genus),
  Outcome_Sex = unique(final_df$Outcome_Sex),
  stringsAsFactors = FALSE
)

heatmap_data <- left_join(all_combinations, final_df %>% select(Genus, Outcome_Sex, Direction),
                          by = c("Genus", "Outcome_Sex")) %>%
  mutate(Direction = ifelse(is.na(Direction), "None", Direction))

heatmap_data$Outcome_Sex <- factor(
  heatmap_data$Outcome_Sex,
  levels = c(
    "FEV1 Z-score (men)", "FEV1 Z-score (women)",
    "FVC Z-score (men)", "FVC Z-score (women)",
    "FEV1/FVC Z-score (men)", "FEV1/FVC Z-score (women)",
    "FeNO (men)", "FeNO (women)"
  )
)

# Plot Figure 3
Figure.3 <- ggplot(heatmap_data, aes(x = Outcome_Sex, y = Genus, fill = Direction)) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Positive" = "#0072B2", "Negative" = "#E69F00", "None" = "white"),
    breaks = c("Positive", "Negative")
  ) +
  labs(x = "", y = "Genus", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(filename = "Figure 3.jpeg", plot = Figure.3, height = 17, width = 14, units = 'in', dpi = 300)
