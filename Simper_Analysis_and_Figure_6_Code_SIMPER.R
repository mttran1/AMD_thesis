#set working directory
setwd("/Users/michelletran/Documents/AMD_thesis")

#load packages
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(vegan)
library(ggplot2)
library(tibble)
library(stats)
library(DirichletMultinomial)
library(DT)
library(stringr)
library(readr)


# Load relative abundance file to run SIMPER ------------------------------
#load otu count matrix # strip "g__" prefixes
otu_counts_og_new_mod <-read_tsv("~/Downloads/FINAL_RERUN_052725/R_folderforanalysis/new_results/otu_table_ra.tsv") %>%
  column_to_rownames(var = "sample_id")
colnames(otu_counts_og_new_mod) <- str_remove(colnames(otu_counts_og_new_mod), "^g__")

otu_counts_og_new_mod <- as.matrix(otu_counts_og_new_mod)  #base matrix for merging with community type data for simper analysis


# read in DMN fit file and get community type assignments -----------------
#load fit file with DMN results and then determine best fit
fit <- readRDS("~/Downloads/FINAL_RERUN_052725/R_folderforanalysis/new_results/re-runDMN/fit_relabund_results_count.rds")
lplc <- sapply(fit, laplace)
plot(lplc, type= "b")
fit[[which.min(lplc)]]
#picking optimal fit
best <- fit[[which.min(lplc)]]
mixturewt(best)

#checking model fit with different number of mixture components using standard information criteria
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
#plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit", )
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

# Plot the lplc line with type "b" (points and lines)

plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit", col = "black", lwd = 2)

# Add the aic line with a different line type
lines(aic, type = "b", lty = 2, col = "blue", lwd = 2)

# Add the bic line with another line type
lines(bic, type = "b", lty = 3, col = "red", lwd = 2)

# Add a legend to describe the three lines
legend(
  "topright",                      # Position of the legend
  legend = c("LPLC", "AIC", "BIC"), # Labels for each line
  col = c("black", "blue", "red"),  # Line colors
  lty = c(1, 2, 3),                 # Line types corresponding to each line
  lwd = 2                           # Line width for the legend59
)


#create a matrix for running DMN data
count <- as.matrix(otu_counts_og_new_mod)

#Extract discrete community assignment for each sample
#    (1-4) from the best-fit DMN
#    (126 samples × K clusters)
post_mat <- fit[[4]]@group
# Confirm rows align with samples
stopifnot(nrow(post_mat) == nrow(count))
# Name rows for clarity
rownames(post_mat) <- rownames(count)

# Get discrete community assignments by highest posterior for each sample
comm_assign <- max.col(post_mat, ties.method = "first")
names(comm_assign) <- rownames(post_mat)
# Check distribution
print(table(comm_assign))

rownames(count) <- rownames(otu_counts_og_new_mod)


# Reorder OTU matrix to match assignment vector
otu2 <- count[names(comm_assign), , drop = FALSE]
stopifnot(all(rownames(otu2) == names(comm_assign)))

# Convert assignments into a factor with informative labels
group_factor <- factor(
  comm_assign,
  levels = c(1, 2, 3, 4),
  labels = c("Community_type_1", "Community_type_2", "Community_type_3", "Community_type_4")
)

levels(group_factor)
print(group_factor)

write_rds(otu2, "otu2.rds")

#calculating average dissimilarity for each pairwise combination using meandist in vegan
meandist(dist = vegdist(otu2), grouping = group_factor)



# Run SIMPER analysis (Bray–Curtis dissimilarity, 999 permutations) -------
sim_out <- simper(
  otu2,
  group_factor,
  permutations = 999
)

# Inspect the top 10 OTU contributors to dissimilarity
sim_sum <- summary(sim_out, ordered = TRUE)
print(sim_sum[[1]][1:10, ])

#Save SIMPER Results for all data
sim_df <- as.data.frame(sim_out[[1]]) %>%
  tibble::rownames_to_column("OTU") %>%
  arrange(desc(average)) %>%
  mutate(cumulative = cumsum(average)) %>%
  select(-species) %>%    # drop the unneeded duplicate
  dplyr::rename(Genus = OTU)

write_csv(sim_df, "SIMPER_Community_types_all.csv")


# Organizing results by comparisons between community types ---------------
# 1) find the names of your comparisons in summary(sim_out)
#    usually there’s just one, like "Community_type_1_Community_type_2"
comparisons <- names( summary(sim_out) )
# e.g. comparisons = "Community_type_1_Community_type_2"

# 2) initialize an empty data.frame
simper_results <- tibble()

# 3) loop (or use purrr::map_dfr) to build the combined table
for (i in 1:length(comparisons)) {
  
  # pull out that comparison and coerce to df
  temp <- summary(sim_out)[[i]] %>% 
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons[i], ".", sep = ""), "", colnames(temp))
  
  # add Comparison & Position & OTU columns
  temp <- temp %>% 
    mutate(
      Comparison = comparisons[i],
      Position   = row_number()) %>%
    rownames_to_column("OTU") 
  # stack
  simper_results <- bind_rows(simper_results, temp) 
}

# 4) inspect
simper_results <- simper_results %>%
  dplyr::rename(Genus = OTU)   # give the same “Genus” name

#5) save comparison results
write_csv(simper_results, "SIMPER_new_ct_comparisons.csv")


# all statistically significant simper results
stat_sig_simper_all <- simper_results %>%
  filter(p <= 0.05)
write_excel_csv(stat_sig_simper_all, "all_stat_sig_results_simper.csv")

# all statistically significant simper results filtered for specific columns
stat_sig_simper_results <- simper_results %>%
  filter(p <= 0.05) %>%
  select(Genus, average, cumsum, Comparison, p, Position) 
write_excel_csv(stat_sig_simper_results, "Final_stat_sig_simper.csv")

# sum of comparisons for all genera
comparison_simper_results <- simper_results %>%
  group_by(Comparison) %>%
  summarize(sum.average = sum(average))

##END OF SIMPER ANALYSIS###


# Sorting SIMPER Results for by community type ----------------------------
library(readr)
library(dplyr)
library(tidyr)

df5 <- read_csv("SIMPER_new_ct_comparisons.csv")


df6 <- df5 %>%
  # 1) extract the two community‐type numbers into new columns CTa and CTb
  extract(
    Comparison,
    into      = c("CTa","CTb"),
    regex     = "Community_type_(\\d+)_Community_type_(\\d+)",
    convert   = TRUE
  ) %>%
  # 2) make a new ‘presence_in_CT’ column
  mutate(
    presence_in_CT = case_when(
      # present only in the first CT (ava ≥ threshold AND avb < threshold)
      ava >= 0.0001 & avb < 0.0001 ~
        paste0("Present in CT", CTa, ", absent in CT", CTb),
      
      # present only in the second CT
      avb >= 0.0001 & ava < 0.0001 ~
        paste0("Present in CT", CTb, ", absent in CT", CTa),
      
      # present in both
      ava >= 0.0001 & avb >= 0.0001 ~
        "Present in both CTs",
      
      # (optional) absent in both
      TRUE ~
        "Absent in both CTs"
    ),
    # 3) —and (optionally) another column for “which is more abundant”
    more_abundant_in = case_when(
      ava >  avb ~ paste0("CT", CTa),
      avb >  ava ~ paste0("CT", CTb),
      TRUE       ~ "equal"
    )
  )

write_csv(df6, "Presence_SIMPER_all_results.csv") #file we use to determine relative abundance for each genus in next step for stat. sig. tests


# filtering through simper results for statistically significant data ----------
significant_simper_5_percent <- tibble(df6) %>%
  filter(p <= 0.05)
  
write_csv(significant_simper_5_percent, "Stat_Sig_simper_relabund.csv")

# Looking at statistically significant simper results ---------------------
#created 5 percent file with relative abundance of genera for top 20% of genera contributing to dissimilarity 
# 1) Read in your full sample × genus abundance table
otu_full <- readRDS("otu2.rds")

# If it's not already a data.frame, coerce it:
otu_full <- as.data.frame(otu_full)

# Then lift the rownames into a real column
otu_full <- otu_full %>%
  rownames_to_column(var = "sample_id")

# Check
head(otu_full)

# 2) Read in your five_percent table and grab the genus list
five_percent <- read_excel(
  "~/Downloads/FINAL_RERUN_052725/R_folderforanalysis/new_results/re-runDMN/5_percent_simper.xlsx"
)

genera_of_interest <- five_percent$Genus

# Subset to genera of interest (vector `genera_of_interest`)
otu_sub <- otu_full %>%
  select(sample_id, all_of(genera_of_interest))

# Load a metadata file that maps each sample to a community type
#    (you’ll need a two‐column file: sample_id, CommunityType with values CT1–CT4)
fit <- readRDS("~/Downloads/FINAL_RERUN_052725/R_folderforanalysis/new_results/re-runDMN/fit_relabund_results_count.rds")
group_assignments_discrete <- factor(apply(fit[[4]]@group, 1, which.max))

# --- Load ecological metadata ---
ecological_data <- read_xlsx("~/Desktop/R_folderforanalysis/new_results/ecological_data_dmn.xlsx") %>%
  column_to_rownames(var = "sample_id")

ecological_data$CT <- group_assignments_discrete

meta <- ecological_data %>%
  # 1) Move rownames into a column called "sample_id"
  rownames_to_column(var = "sample_id") %>%
  # 2) Keep only sample_id + CT
  select(sample_id, CT)


# 5) Melt to long form and join metadata
df_long <- otu_sub %>%
  pivot_longer(
    cols      = -sample_id,
    names_to  = "Genus",
    values_to = "Abundance"
  ) %>%
  left_join(meta, by = "sample_id") 
# now df_long has columns: sample_id, Genus, Abundance, CT

# Run per‐genus ANOVA & Kruskal–Wallis
#    and Bonferroni‐adjust the p‐values
anova_res <- df_long %>%
  group_by(Genus) %>%
  summarise(
    p_anova = summary(aov(Abundance ~ CT, data = .))[[1]][["Pr(>F)"]][1]
  ) %>%
  ungroup() %>%
  mutate(p_anova_bonf = p.adjust(p_anova, method = "bonferroni"))

kw_res <- df_long %>%
  group_by(Genus) %>%
  summarise(
    p_kw = kruskal.test(Abundance ~ CT, data = .)$p.value
  ) %>%
  ungroup() %>%
  mutate(p_kw_bonf = p.adjust(p_kw, method = "bonferroni"))

# Show significant genera
sig_anova   <- filter(anova_res, p_anova_bonf < 0.05)
sig_kruskal <- filter(kw_res,   p_kw_bonf   < 0.05)

print(sig_anova)
print(sig_kruskal)


#(Optional) Dunn’s test for pairwise comparisons on those significant genera
#    install.packages("FSA") if you haven’t already
library(FSA)
dunn_list <- df_long %>%
  filter(Genus %in% sig_kruskal$Genus) %>%
  group_by(Genus) %>%
  do(dunn = dunnTest(Abundance ~ CT, data = ., method="bonferroni")$res) %>%
  unnest(cols = c(dunn))

print(dunn_list)


#making df_long with CT and genera of interest
#Filter & summarize
library(dplyr)
library(tidyr)

summary_df <- df_long %>%
  filter(Genus %in% genera_of_interest) %>%
  group_by(Genus, CT) %>%
  summarise(relative_abundance = sum(Abundance, na.rm = TRUE), .groups="drop")

write_csv(summary_df, "summary_5percent_ct_abund.csv")

# Pivot to wide format so each CT is its own column
wide_df <- summary_df %>%
  pivot_wider(
    names_from   = CT,
    values_from  = relative_abundance,
    names_prefix = "CT",
    values_fill  = 0
  ) %>%
  arrange(Genus)

# View your table
print(wide_df)
write_csv(wide_df, "5percent_by_ct.csv")

# running anova and tests for sig of 5percent -----------------------------

# 1) Filter your sample‐level long data
df_sub <- df_long %>%
  filter(Genus %in% genera_of_interest) %>%
  mutate(CT = factor(CT))    # ensure CT is factor

# 2) ANOVA per genus
anova_sub <- df_sub %>%
  group_by(Genus) %>%
  summarise(
    p_anova      = summary(aov(Abundance ~ CT, data = .))[[1]][["Pr(>F)"]][1]
  ) %>%
  ungroup() %>%
  mutate(p_anova_bonf = p.adjust(p_anova, method = "bonferroni"))

# 3) Kruskal–Wallis per genus
kw_sub <- df_sub %>%
  group_by(Genus) %>%
  summarise(
    p_kw      = kruskal.test(Abundance ~ CT, data = .)$p.value
  ) %>%
  ungroup() %>%
  mutate(p_kw_bonf = p.adjust(p_kw, method = "bonferroni"))

# 4) Merge back into your summary_df (which has Genus, CT, total_abundance)
summary_with_p <- summary_df %>%
  filter(Genus %in% genera_of_interest) %>%      # keep only those rows
  left_join(kw_sub,    by = "Genus") %>%
  filter(p_kw_bonf < 0.05)             # keep only those with adjusted p < 0.05

# 5) Inspect
print(summary_with_p)

write_csv(summary_with_p, "analysis_sig_no_anova_5percent.csv")

m <- aov(Abundance ~ CT, data = df_sub)  
shapiro.test(residuals(m))
qqnorm(residuals(m)); qqline(residuals(m))
#normality tests all violated so no doing anova

kruskal.test(Abundance ~ CT, data = df_sub)
#Extremely significant difference in median abundance across your four CTs for this genus. You reject the null hypothesis that all CT distributions are the same.
library(car)
leveneTest(Abundance ~ CT, data = df_sub)
#No evidence of unequal variances (p > 0.05), so the homogeneity‐of‐variance assumption is met.

#getting adjusted p-values 
library(FSA)  
dunn_res <- dunnTest(Abundance ~ CT, data = df_sub, method = "bonferroni")  

#testing effect size
library(rstatix)
df_sub %>% kruskal_effsize(Abundance ~ CT)

print(dunn_res)

# 1) Filter to your genera of interest
df_sub <- df_long %>% 
  filter(Genus %in% genera_of_interest)

# 1) Omnibus Kruskal–Wallis + Bonferroni + stars
kw_summary <- df_sub %>%
  group_by(Genus) %>%
  kruskal_test(Abundance ~ CT) %>%          # χ², df, p
  adjust_pvalue(method = "bonferroni") %>%  # p.adj
  add_significance("p.adj")                 # p.adj.signif (“ns”, “*”, etc.)

# 2) Effect sizes η²[H]
effsize_df <- df_sub %>%
  group_by(Genus) %>%
  kruskal_effsize(Abundance ~ CT)           # gives effsize & magnitude

# 3) Join them
kw_summary <- left_join(kw_summary, effsize_df, by = "Genus")

# 4) Inspect
kw_summary

#cleaning up because of duplicate info
kw_clean <- kw_summary %>%
  select(
    Genus,
    statistic,        # χ²
    df,               # degrees of freedom
    p,                # raw p
    p.adj,            # Bonferroni‐adjusted p
    p.adj.signif,     # “ns”, “*”, etc.
    effsize,          # η²[H]
    magnitude         # small/medium/large
  )

# Inspect
print(kw_clean)
write_csv(kw_clean, "kw_clean_5_percent_final.csv")

#kw test adj. p-values failed 
kw_sub <- kw_sub %>%
  mutate(p_kw_fdr = p.adjust(p_kw, method = "BH"))  # Benjamini-Hochberg FDR correction

#joining them together
df_sub_with_p <- df_sub %>%
  left_join(
    kw_sub %>% select(Genus, p_kw, p_kw_fdr), 
    by = "Genus"
  )

significant_df <- df_sub_with_p %>%
  filter(p_kw_fdr < 0.05)

# 1. Calculate adjusted p-values on kw_sub
kw_sub <- kw_sub %>%
  mutate(p_kw_fdr = p.adjust(p_kw, method = "BH"))

# 2. Join to df_sub by Genus
df_sub_with_p <- df_sub %>%
  left_join(kw_sub %>% select(Genus, p_kw, p_kw_fdr), by = "Genus")

# 3. Filter or plot using adjusted p-values
significant_df <- df_sub_with_p %>%
  filter(p_kw_fdr < 0.05)

write_csv(significant_df, "analysis_sig_kw_fdr_5percent.csv")
# creating Figure 6 for Genus relative abundance v. ct box plot -----------
library(tidyverse)
library(ggpubr)
library(readr)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)

fivepercent_analysis_summ <- read.csv("analysis_sig_kw_fdr_5percent.csv")
wide_abundances <- fivepercent_analysis_summ 
    
    
# Step 1: Load your original p-value table
pval_df <- read_csv("Presence_SIMPER_all_results.csv")
# Step 2: Create a column name for each comparison
pval_df <- pval_df %>%
  mutate(p_col = paste0("p_", CTa, CTb))

# Step 3: before pivot remove duplicate p-values
pval_unique <- pval_df %>%
  group_by(Genus, p_col) %>%
  summarise(p = dplyr::first(p), .groups = "drop")

# Step 4: Pivot wider to get each p_ij in its own column
pval_wide <- pval_unique %>%
  select(Genus, p_col, p) %>%
  pivot_wider(names_from = p_col, values_from = p)

# Step 5: Join with wide_abundances
wide_abundances_with_p <- wide_abundances %>%
  left_join(pval_wide, by = "Genus")

write_csv(wide_abundances_with_p, "CT_all_abund_for_figure.csv")

# 0) Your palette (names must exactly match recoded groups)
my_cols <- c(
  "Community type 1" = "#F8766D",
  "Community type 2" = "#7CAE00",
  "Community type 3" = "#00BFC4",
  "Community type 4" = "#C77CFF"
)


df_long2 <- read_csv("CT_all_abund_for_figure.csv") 
df_long2 <- df_long2 %>%
  select(-p_kw, -p_kw_fdr)#removing p_kw and p_kw_fdr columns


# Step 1: Ensure df_long2$CT is a factor with the correct level order
df_long2 <- df_long2 %>%
  mutate(CT = factor(CT, levels = 1:4, labels = paste("Community type", 1:4)))

pval_long <- df_long2 %>%
  distinct(Genus, .keep_all = TRUE) %>%
  select(Genus, starts_with("p_")) %>%
  pivot_longer(cols = starts_with("p_"),
               names_to = "comparison",
               values_to = "p_val") %>%
  mutate(
    CT1 = paste0("Community type ", substr(comparison, 3, 3)),
    CT2 = paste0("Community type ", substr(comparison, 4, 4)),
    sig_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

#adding significance bars
pval_long <- pval_long %>%
  mutate(
    CT1 = str_trim(CT1),
    CT2 = str_trim(CT2),
    CT1_num = as.numeric(factor(CT1, levels = paste("Community type", 1:4))),
    CT2_num = as.numeric(factor(CT2, levels = paste("Community type", 1:4))),
    x       = (CT1_num + CT2_num) / 2
  )

max_heights <- df_long2 %>%
  group_by(Genus) %>%
  summarise(max_y = max(Abundance, na.rm = TRUE), .groups = "drop")

pval_long <- pval_long %>%
  left_join(max_heights, by = "Genus") %>%
  group_by(Genus) %>%
  arrange(CT1_num, CT2_num) %>%
  mutate(
    y_bracket = max_y + 0.05 * max_y * row_number(),  # bracket height
    y_label   = y_bracket + 0.02 * max_y              # text height
  ) %>%
  ungroup()



# Plot
ggplot(df_long2, aes(x = CT, y = Abundance, fill = CT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  # Add significance brackets (lines)
  geom_segment(
    data = pval_long %>% filter(sig_label != "ns"),
    aes(
      x = CT1_num,
      xend = CT2_num,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  # Add significance stars or labels
  geom_text(
    data = pval_long %>% filter(sig_label != "ns"),
    aes(
      x = x,
      y = y_label,
      label = sig_label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10)
  ) +
  labs(
    x = "Community Type",
    y = "Relative Abundance",
    fill = "Community Type",
    caption = "*** p < 0.001   ** p < 0.01   * p < 0.05"
  )




