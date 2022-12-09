# Preps ####

# Load library
library(BiodiversityR)
library(SpadeR)
library(iNEXT)
library(treemap)
library(cluster)
library(tidyverse) 
library(abdiv)

# Load dataset


# Make tabel of diversity index ####

# Reorder dataframe for calculating species diversity
mi.spabd <- raw %>% count(Scientific.name, Transect) %>% 
  pivot_wider(names_from = Scientific.name, values_from = n)
mi.spabd[is.na(mi.spabd)] <- 0
mi.df <- column_to_rownames(mi.spabd, var="Transect")

# Create a list of index from BiodiversityR packages
index_list <- c("richness", "abundance", "Shannon", "Jevenness")

# Run a loop over the list of indices
result <- lapply(index_list, function(x) {
  # Calculate the diversity index for each index
  result <- diversityresult(mi.df, index = x, method = "each site", digits = 2)
  
  # Return the result
  return(result)
})

# Print the results
print(result)

# Merge the results
result_d.i <- do.call(cbind, result)

# add margalef and dominance index from abdiv package
Margalef <- mi.df %>% t() %>% apply(2, margalef)
Dominance <- mi.df %>% t() %>% apply(2, dominance)

#merge all index
(indeks_tabel <- cbind(result_df, Margalef, Dominance))

## Make plots from diversity index ####

# richness plot
richness_plot <- indeks_tabel %>% select((richness:abundance)) %>% rownames_to_column(var = "Transect") %>%
  pivot_longer(!Transect, names_to = "index", values_to = "values")

ggplot(richness_plot, aes(fill=index, y=values, x=Transect,)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label = round(values, 1)), 
            position = position_dodge(0.9), vjust = -0.5, hjust = 0.5)


# index plot
index_plot <- indeks_tabel %>% select(!(richness:abundance)) %>% rownames_to_column(var = "Transect") %>%
  pivot_longer(!Transect, names_to = "index", values_to = "values")

ggplot(index_plot, aes(fill=Transect, y=values, x=index,)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label = round(values, 1)), 
            position = position_dodge(0.9), vjust = -0.5, hjust = 0.5)


# Rarefaction and estimating species richnes ####
mi.abd <- raw %>% count(Scientific.Name) %>%  select (n)

# estimating species richness
out1 <- SpadeR::ChaoSpecies(mi.abd, datatype = "abundance", k=10,conf=0.95)
out1$Basic_data_information
out1$Species_table

# produce plot
out2 <- iNEXT::iNEXT(mi.abd, q=0, datatype="abundance")
ggiNEXT(x=out2, type=1) #Sample-size-based rarefaction
ggiNEXT(x=out2, type=2) #Sample completeness curve
ggiNEXT(x=out2, type=3) #Coverage-based rarefaction


# Clustering sample or transect #### 

# calculate dissimilarity from jaccard 
distmatrix <- vegdist(mi.df, method="jaccard")
Clust.1 <- diana(distmatrix)
summary(Clust.1)

#produce plot
plot(Clust.1, which.plots=2,hang=-1,main="",sub="",xlab="",ylab="")
