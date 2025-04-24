#set working directory 

setwd("~/Desktop/R_folderforanalysis/new_results")

# workflow taken from https://bioconductor.org/packages/release/bi --------
library("DT")
library(DirichletMultinomial)
library(lattice)
library(parallel)
library(readr)
library(readxl)
library(dplyr)
library(phyloseq)
library(microbiome)
library(reshape2)
library(magrittr)
library(tidyr)

#load fit file
fit <- readRDS("~/Desktop/R_folderforanalysis/new_results/fit_results_count.rds")

lplc <- sapply(fit, laplace)
plot(lplc, type= "b")
fit[[which.min(lplc)]]

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
  lwd = 2                           # Line width for the legend
)

#picking optimal fit
best <- fit[[which.min(lplc)]]
mixturewt(best)

for (k in 1:2) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  
  d2 <- d %>% 
    filter(cluster == k) %>%
    mutate(OTU = gsub("^g__", "", OTU)) %>%             # strip g__ prefix
    arrange(desc(abs(value))) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    slice_max(order_by = abs(value), n = 30)
  
  p2 <- ggplot(d2, aes(x = OTU, y = value)) +
    geom_col() +
    coord_flip() +
    labs(
      tag   = letters[k],                                # Adds a / b
      x     = "OTU",
      y     = "Value (unitless α)"
    ) +
    theme_minimal() +
    theme(
      plot.tag.position = c(0.02, 0.98),                 # top-left corner
      plot.tag           = element_text(size = 16, face = "bold"),
      axis.text.y        = element_text(size = 8),
      plot.title         = element_text(size = 14, face = "bold"),
      axis.title         = element_text(size = 12)
    )
  
  print(p2)
  ggsave(
    paste0("top_drivers_community_type_", k, ".png"),
    plot   = p2,
    width  = 8, height = 6
  )
}

#taxa fitted to Dirichlet Components
splom(log(fitted(best)))
dev.off()

#posterior mean diff
pO <- fitted(fit[[1]], scale=TRUE) #scale by theta
p2 <- fitted(best, scale=TRUE)
colnames(p2) <- paste("m", 1:2, sep="")
(meandiff <- colSums(abs(p2 - as.vector(pO))))
sum(meandiff)

#to produce taxonomic contribution (10 largest) to Dirichlet components
diff <- rowSums(abs(p2 - as.vector(pO)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o] / sum(diff))
df <- head(cbind(Mean=pO[o], p2[o], diff=diff[o], cdiff), 30)
DT::datatable(df) |>
  DT::formatRound(c("Mean", "V2", "diff", "cdiff"), digits = 4)

