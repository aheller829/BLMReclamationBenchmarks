# Fuzzy cluster analysis example for ecoregion 24
# Load packages
library(vegan)
library(cluster)
library(dplyr)
library(tidyverse)
library(NbClust)
library(RColorBrewer)

# Read in monitoring data tables
l3.indicators <- read.csv("Example\\l3.24.indicators.csv")
l3.plots <- read.csv("Example\\l3.24.plots.csv")

set.seed(1)

# Define analysis table and add environmental variables to indicators
analysis.table <- l3.indicators
names(l3.plots)
names(l3.indicators)
env <- dplyr::select(l3.plots, PrimaryKey, ElevationM, PRISM)
env <- subset(env, env$PrimaryKey %in% analysis.table$PrimaryKey)
analysis.table <- dplyr::left_join(analysis.table, env)
# Remove NAs
analysis.table <- na.omit(analysis.table)
env <- subset(env, env$PrimaryKey %in% analysis.table$PrimaryKey)
# Add numeric id
analysis.table$ID <- rownames(analysis.table)
analysis.table$ID <- as.numeric(analysis.table$ID)


# Build dissimilarity matrix using indicators
names(analysis.table)
# Select indicators
fg.foliar <- analysis.table %>%
  dplyr::select(S, GapLg, TotalFoliar, Nox, Inv, max_height, NatHerb, NatWoody) 
# Convert to a gower dissimilarity matrix
fg.foliar.dist <- vegan::vegdist(fg.foliar, method = "gower", binary = FALSE)

# Run a PCoA on functional groups and structural indicators
vegPCA <- cmdscale(fg.foliar.dist, k = 2)

# View ordination
ordiplot(vegPCA, type = "text")
plot(vegPCA)
# Look at functional group correlation with ordination axes and plot over ordination
fit <- envfit(vegPCA, fg.foliar, perm = 999, na.rm = TRUE, choices = c(1, 2, 3))
fit
plot(fit, p.max = 0.05, col = "red")
# Look at environmental variable correlation with ordination axes
env <- dplyr::select(env, -PrimaryKey)
plot(vegPCA)
fit <- vegan::envfit(vegPCA, env, perm = 999, na.rm = TRUE, choices = c(1, 2, 3))
fit
plot(fit, p.max = 0.05, col = "blue")



# Fuzzy clustering
# Start by finding optimal number of clusters indicated by cluster metrics (k)

# Test cluster number metrics
fg.KM.cascade <- vegan::cascadeKM(fg.foliar.dist, inf.gr = 2, sup.gr = 15, iter = 100, criterion = "ssi")
plot(fg.KM.cascade, sortg = TRUE)

# Adjust fuzziness/crispness with membership exponent approaching 2 for fuzzier classification
veg.fanny <- cluster::fanny(fg.foliar.dist, k = 11, memb.exp = 1.1, maxit = 1000, keep.diss = TRUE)
# Display's Dunn's partition coefficient (low coeff = very fuzzy, near 1 = crisp)
veg.fanny$coeff
# Build a dataframe of membership values
fanny.mems <- as.data.frame(veg.fanny$membership)
fanny.mems <- fanny.mems %>%
  dplyr::mutate_if(is.numeric, round, digits = 3)

# Plot fuzzy clusters in ordination plot
# Set color palette
colors <- brewer.pal(n = 12,"Set3")
plot(vegPCA)
stars(veg.fanny$membership, locatio = vegPCA, draw.segm = TRUE, add = TRUE, scale = FALSE, len = 0.05, 
      col.segments = colors, labels = NULL)
ordihull(vegPCA, veg.fanny$clustering, col = "black")
ordispider(vegPCA, veg.fanny$clustering, col = "gray", label = T)


# Assign plots to clusters by top membership value
names(fanny.mems)
topmems <- fanny.mems %>%
  dplyr::mutate(ID = rownames(fanny.mems)) %>%
  tidyr::gather(Cluster, MemVal, V1:V11) %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(MemVal) %>%
  dplyr::slice(which.max(MemVal)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Cluster, MemVal)
str(analysis.table)
str(topmems)
topmems$ID <- as.numeric(topmems$ID)
topmems<- inner_join(topmems, analysis.table, by = "ID")

# Summarize indicator means across clusters
names(topmems)
topmem.summary <- topmems %>%
  dplyr::select(Annual_Forb:max_height, Cluster) %>%
  tidyr::gather(Species, MeanCover, Annual_Forb:max_height, factor_key = FALSE) %>%
  dplyr::group_by(Cluster, Species) %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::filter(MeanCover > 0)
# Sort by highest
# Turn off scientific notation
options(scipen = 999)
topmem.summary <- topmem.summary %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(desc(MeanCover), .by_group = TRUE)

# Separate plots that have high membership values to a cluster from "fuzzy" plots
high.mems <- filter(topmems, MemVal > 0.50)
low.mems <- filter(topmems, MemVal < 0.51)


# Save to csv
write.csv(high.mems, "high.mems.er.24.csv", row.names = FALSE)



