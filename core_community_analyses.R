
#identifying core community from sourmash results using relative abundance--------


#setting working directory
setwd("~/Desktop/R_folderforanalysis/new_results")

# install.packages(c("readxl","dplyr","purrr"))
library(readxl)
library(dplyr)
library(purrr)


# Checking core community at 95% presence across samples ------------------

# 1. Compute your df_95 as before:
df <- read_xlsx("otu_table_ra.xlsx")
presence <- df %>% mutate(across(-1, ~ .x > 0))
n_samples <- nrow(df)
genus_95 <- presence %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x))) %>%
  pivot_longer(everything(), names_to="genus", values_to="n_present") %>%
  filter(n_present >= 0.95 * n_samples) %>%
  pull(genus)

df_95 <- df %>% select(1, all_of(genus_95))

# 2. Summarise counts for those 95% genera:
counts_95 <- df_95 %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# 3. Build your “count row” tibble:
count_row_95 <- tibble(run_accessions = "n_genera_present") %>%
  bind_cols(counts_95)

# 4. Append it:
df_95_with_count <- bind_rows(df_95, count_row_95)


# Checking core community at 90% presence across samples --------------


# 1. Compute your df_90 as before:
df <- read_xlsx("otu_table_ra.xlsx")
presence <- df %>% mutate(across(-1, ~ .x > 0))
n_samples <- nrow(df)
genus_90 <- presence %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x))) %>%
  pivot_longer(everything(), names_to="genus", values_to="n_present") %>%
  filter(n_present >= 0.90 * n_samples) %>%
  pull(genus)

df_90 <- df %>% select(1, all_of(genus_90))

# 2. Summarise counts for those 90% genera:
counts_90 <- df_90 %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# 3. Build your “count row” tibble:
count_row_90 <- tibble(run_accessions = "n_genera_present") %>%
  bind_cols(counts_90)

# 4. Append it:
df_90_with_count <- bind_rows(df_90, count_row_90)


# Checking core community at 98% presence across samples ------------------


# 1. Compute your df_98 as before:
df <- read_xlsx("otu_table_ra.xlsx")
presence <- df %>% mutate(across(-1, ~ .x > 0))
n_samples <- nrow(df)
genus_98 <- presence %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x))) %>%
  pivot_longer(everything(), names_to="genus", values_to="n_present") %>%
  filter(n_present >= 0.98 * n_samples) %>%
  pull(genus)

df_98 <- df %>% select(1, all_of(genus_98))

# 2. Summarise counts for those 98% genera:
counts_98 <- df_98 %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# 3. Build your “count row” tibble:
count_row_98 <- tibble(run_accessions = "n_genera_present") %>%
  bind_cols(counts_98)

# 4. Append it:
df_98_with_count <- bind_rows(df_98, count_row_98)




# identifying the core communit from sourmash otu counts-> 95%--------


# 1. Compute your df_95 as before:
df_otu <- read_xlsx("~/Desktop/R_folderforanalysis/OTU_table_genus_final.xlsx")
presence <- df_otu %>% mutate(across(-1, ~ .x > 0))
n_samples <- nrow(df_otu)
genus_otu_95 <- presence %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x))) %>%
  pivot_longer(everything(), names_to="genus", values_to="n_present") %>%
  filter(n_present >= 0.95 * n_samples) %>%
  pull(genus)

df_otu_95 <- df_otu %>% select(1, all_of(genus_otu_95))

# 2. Summarise counts for those 95% genera:
counts_otu_95 <- df_otu_95 %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# 3. Build your “count row” tibble:
count_otu_row_95 <- tibble(run_accessions = "n_genera_present") %>%
  bind_cols(counts_otu_95)

# 4. Append it:
df_otu_95_with_count <- bind_rows(df_otu_95, count_otu_row_95)

# identifying the core community from otu counts 90% ----------------------

# 1. Compute your df_95 as before:
df_otu <- read_xlsx("~/Desktop/R_folderforanalysis/OTU_table_genus_final.xlsx")
presence <- df_otu %>% mutate(across(-1, ~ .x > 0))
n_samples <- nrow(df_otu)
genus_otu_90 <- presence %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x))) %>%
  pivot_longer(everything(), names_to="genus", values_to="n_present") %>%
  filter(n_present >= 0.90 * n_samples) %>%
  pull(genus)

df_otu_90 <- df_otu %>% select(1, all_of(genus_otu_90))

# 2. Summarise counts for those 90% genera:
counts_otu_90 <- df_otu_90 %>%
  select(-1) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# 3. Build your “count row” tibble:
count_otu_row_90 <- tibble(run_accessions = "n_genera_present") %>%
  bind_cols(counts_otu_90)

# 4. Append it:
df_otu_90_with_count <- bind_rows(df_otu_90, count_otu_row_90)

#Results show that using OTU counts or sourmash tax metagenome abundance results produce the same microbes identified for 90%, 95%, and 98% categories
# First three code sections produce the number of samples the genera identified can be found across from the 126 samples
