
# Relative abundance bar plot ---------------------------------------------
setwd("~/Desktop/R_folderforanalysis/new_results")
# Required libraries
install.packages("tidytext")
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(gridExtra) #for multiple plots
library(forcats)
library(tidytext)
library(stringr)

# Ensure the environmental data is properly structured
# select best fit (best Large DMN) fit to run NMDS of results with metadata: pH, Fe[II]; mine type; country--------
metadata <- read_csv("~/Desktop/R_folderforanalysis/new_results/cleaned_modified_full_metadata.csv")

# Remove empty rows (handles different column types)
metadata_clean <- metadata[apply(metadata, 1, function(row) {
  all(is.na(row) | row == "" | trimws(as.character(row)) == "")
}) == FALSE, ]

# View cleaned data
head(metadata_clean)

otu_counts_og_new_mod <- read_excel("~/Desktop/R_folderforanalysis/OTU_table_genus_final.xlsx")

# pH nmds of DMN data -----------------------------------------------------
#add other environmental variables, [Fe], ORP, mine type, country
metadata_env_var <- metadata_clean %>%
  select('run_accessions', 'pH', 'isolation_source', 'mine_type', 'Geographic_location', '[Fe(II)]_(mg/kg)', '[Fe(II)]_(mg/L)_in_water', 'Total_C_ TC_(mg/L)_in_water', '[SO42-]_(mg/kg)')

# Clean run_accessions columns in both data frames
otu_counts_og_new_mod$run_accessions <- tolower(trimws(as.character(otu_counts_og_new_mod$run_accessions)))
metadata_env_var$run_accessions <- tolower(trimws(as.character(metadata_env_var$run_accessions)))

#binding env_var table with otu matrix
env_otu_matrix <- left_join(metadata_env_var, otu_counts_og_new_mod, by = "run_accessions")


#write env_otu table to folder
envcounts_matrix <- file.path("~/Desktop/R_folderforanalysis", "env_counts_matrix.tsv") 
write.table(env_otu_matrix, file = envcounts_matrix, sep = "\t", row.names = FALSE)

env_data <- read_tsv("~/Desktop/R_folderforanalysis/new_results/env_otu_matrix.tsv")
env_data <- env_data %>%
  rename_with(
    ~ str_remove(., "^g__"),    # function to apply to names
    starts_with("g__")          # only on columns beginning with “g__”
  )


# Add pH information and taxa of interest
env_data$pH <- metadata_env_var$pH  # Assuming pH is a column in your metadata

# Select the genera of interest
genera_of_interest <- c("Ferrovum", "Leptospirillum", "Acidithiobacillus", "Gallionella")

# Filter and summarize data
library(dplyr)
library(ggplot2)
library(forcats)
# Assuming each genus is a column in the fitted matrix
# Step 1: Calculate the sample counts for each pH group
sample_counts <- env_data %>%
  mutate(pH_group = cut(pH, breaks = c(-Inf, 2, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, Inf), 
                        labels = c("<2.0", "2.0-2.4", "2.4-2.8", "2.8-3.2", 
                                   "3.2-3.6", "3.6-4.0", "4.0-4.4", "4.4-4.8", ">4.8"))) %>%
  count(pH_group)  # Count the number of samples in each group

# Step 2: Add sample counts to the x-axis labels
sample_counts <- sample_counts %>%
  mutate(pH_group_label = paste0(pH_group, "\n(n=", n, ")"))  # Create labels with counts

# Step 3: Prepare relative abundances data
relative_abundances <- env_data %>%
  select(pH, all_of(genera_of_interest)) %>%  # Select relevant columns
  pivot_longer(cols = all_of(genera_of_interest), 
               names_to = "Genus", values_to = "Abundance") %>%  # Reshape data
  mutate(pH_group = cut(pH, breaks = c(-Inf, 2, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, Inf), 
                        labels = c("<2.0", "2.0-2.4", "2.4-2.8", "2.8-3.2", 
                                   "3.2-3.6", "3.6-4.0", "4.0-4.4", "4.4-4.8", ">4.8"))) %>%
  group_by(pH_group, Genus) %>%
  summarise(Relative_Abundance = mean(Abundance, na.rm = TRUE) * 100, .groups = "drop")

# Step 4: Merge sample count labels with relative abundances
relative_abundances <- relative_abundances %>%
  left_join(sample_counts, by = "pH_group") %>%
  mutate(pH_group_label = paste0(pH_group, "\n(n=", n, ")"))  # Add sample count labels

# Step 5: Create the bar plot with updated x-axis labels
ggplot(relative_abundances, aes(x = pH_group_label, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "pH (with Sample Counts)", 
    y = "Relative Abundance (%)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("violet", "red", "purple", "orange")) +  # Customize colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Checking relative abundance values, if values are between 0-1 then the rel. abundance is 0-1 but if they are greater it is a percentage
env_data %>%
  select(all_of(genera_of_interest)) %>%
  summarise(across(everything(), range, na.rm=TRUE))



