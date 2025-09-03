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







# creating Figure 6 for Genus relative abundance v. ct box plot -----------
library(tidyverse)
library(ggpubr)
library(dplyr)
library(readr)


# 0) Your palette (names must exactly match recoded groups)
my_cols <- c(
  Community_type_1 = "#F8766D",
  Community_type_2 = "#7CAE00",
  Community_type_3 = "#00BFC4",
  Community_type_4 = "#C77CFF"
)

df <- read_csv

