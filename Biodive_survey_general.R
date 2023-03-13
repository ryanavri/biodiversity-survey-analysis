#Preparation----
##Load library----
library(tidyverse)
library(BiodiversityR)
library(SpadeR)
library(iNEXT)
library(ggrepel)
library(ggplot2)
library(ggvegan)

##Load dataset----

#Analysis of richness and diversity----

##Produce abundance, richness and various diversity index----
indeks_tabel <- raw %>%
  count(Scientific.Name, Transect) %>%
  group_by(Transect) %>%
  summarize(richness = n_distinct(Scientific.Name),
            abundance = sum(n),
            shannon = -sum(prop.table(n) * log(prop.table(n))),
            margalef = (n_distinct(Scientific.Name) - 1) / log(sum(n)),
            evenness = (-sum(prop.table(n) * log(prop.table(n))))/log(length(n)),
            simpson = sum(prop.table(n)^2))%>%
  mutate(across(4:last_col(), ~round(., 2)))

#if you use vegan and well known with the formula, you may use this simplified code
#library(vegan)
#indeks_tabel <- raw %>%
#  count(Scientific.Name, Transect) %>%
#  group_by(Transect) %>%
#  summarize(richness = specnumber(n),
#            abundance = sum(n),
#            shannon = diversity(n, index = "shannon"),
#            margalef = (specnumber(n) - 1) / log(sum(n)),
#            evenness = shannon /log(length(n)),
#            simpson = 1 - diversity(n, index = "simpson"))


###Plot richness and abundance plot----
(richness_plot <- indeks_tabel %>% 
   select(c(Transect, richness, abundance)) %>% 
   pivot_longer(-Transect, names_to = "category", values_to = "values") %>%
   ggplot(aes(fill=category, y=values, x=Transect)) + 
   geom_col(position="dodge", width = 0.8) + 
   geom_text(aes(label = round(values,1)), 
             position = position_dodge(0.8), vjust = -0.5, hjust = 0.5))

###Plot diversity index plot----
(RI_plot <- indeks_tabel %>% 
   select(c(Transect,shannon,margalef)) %>% 
   pivot_longer(-Transect, names_to = "index", values_to = "values") %>%
   ggplot(aes(fill=index, y=values, x=Transect)) + 
   geom_col(position="dodge", width = 0.8) + 
   geom_text(aes(label = round(values,2)), 
             position = position_dodge(0.8), vjust = -0.5, hjust = 0.5))

###Plot evenness index plot----
(DI_plot <- indeks_tabel %>% 
   select(c(Transect,simpson,evenness)) %>% 
   pivot_longer(-Transect, names_to = "index", values_to = "values") %>%
   ggplot(aes(fill=index, y=values, x=Transect)) + 
   geom_col(position="dodge", width = 0.8) + 
   geom_text(aes(label = round(values,2)), 
             position = position_dodge(0.8), vjust = -0.5, hjust = 0.5))

##Count of species----

###Plot the abundance of each species----

#produce species abundance to reorder by the highest abd. species
species_abundance <- raw_m %>%
  group_by(Scientific.Name) %>%
  summarise(abundance = n()) %>%
  arrange(desc(abundance))

#plot the abundance of each species
ggplot(raw_m, aes(x = factor(Scientific.Name, levels = rev(species_abundance$Scientific.Name)), fill = Transect)) +
  geom_bar() +
  labs(x = "Species", y = "Number of Individuals", fill = "Transect") +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic")) + 
  coord_flip()

#alternatif fancy plot
ggplot(raw_m, aes(x = factor(Scientific.Name, levels = rev(species_abundance$Scientific.Name)), fill = Transect)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9)) +
  labs(x = "Species", y = "Number of Individuals", fill = "Transect") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic", size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")) +
  guides(fill = guide_legend(title = "Transect", title.position = "top", title.hjust = 0.5)) +
  coord_flip()

###Plot the abundance of dominating species----
abundance_table <- raw_m %>%
  group_by(Transect, Scientific.Name) %>%
  summarise(abundance = n())

ggplot(abundance_table, aes(x = factor(Scientific.Name, levels = rev(sort(unique(Scientific.Name)))), y = abundance, fill = Transect)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Abundance", fill = "Transect") +
  facet_wrap(~ Transect, ncol = 2) +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()

#alternatif fancy plot
ggplot(abundance_table, aes(x = factor(Scientific.Name, levels = rev(sort(unique(Scientific.Name)))), y = abundance, fill = Transect)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Species", y = "Abundance", fill = "Transect") +
  facet_wrap(~ Transect, ncol = 2) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()


###Calculates the importance values of tree species----
IV <- raw_t %>% 
  group_by(Transect, Scientific.Name) %>% 
  summarise(count = n(),
            basal = sum(0.7854*(dbh/100)^2)) %>% #this convert cm dbh into msq for basal area
  ungroup() %>% 
  as.data.frame() %>%
  importancevalue(site="Transect", species="Scientific.Name", 
                  count="count", basal="basal", 
                  factor="",level="") %>%
  as.data.frame() %>%
  mutate(across(3:last_col(), ~round(., 2)))

#Analysis of differences in species composition----

##Analysis of ecological distance by clustering----

# Create a data matrix with species data
data_matrix <-table(raw$Transect, raw$Scientific.Name)

# Calculate the Bray-Curtis dissimilarity matrix and plot the dendogram
data_matrix %>%
  vegdist(method = "bray") %>%
  hclust(method = "average") %>%
  plot(xlab="", ylab="Dissimilarity", sub="Transect", hang=-1)

##Analysis of ecological distance by ordinance----

# Create abundance matrix with species data
sp.abd <- raw %>%
  count(Transect, Scientific.Name) %>%
  pivot_wider(names_from = Scientific.Name, values_from = n, values_fill = list(n = 0)) %>%
  column_to_rownames(var = "Transect")

# Calculate NMDS with bray-curtis
mynmds <- metaMDS(sp.abd, distance = "bray")

# Group the points by nearest distances using k-means clustering
myclusters <- kmeans(mynmds$points, centers=3, nstart=10) #change centers value for different group sizes
group <- factor(myclusters$cluster)

# Plot the NMDS ordination with grouped points and plot names
plot(mynmds, type="n", main="NMDS Plot")
points(mynmds, display="sites", col=group, pch=16)
text(mynmds$points, labels=rownames(sp.abd), cex=0.7, pos = 1)
legend("bottomright", legend=unique(group), col=unique(group), pch=16, title = "Group")

#Interpolation and Extrapolation for Species Diversity----

##Perform the species accumulation curve----
out1 <- raw %>% count(Scientific.Name) %>% 
  column_to_rownames(var = "Scientific.Name") %>% 
  as.data.frame() %>%
  iNEXT(q=0, datatype = "abundance", conf = 0.95, nboot = 100) #the bootstrap can be set to higher number

#plot the result
ggiNEXT(x = out1, type = 1, color.var = "Order.q") +
  labs(title = "Species Accumulation Curve", x = "Number of Individuals", y = "Cumulative Species Richness") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

##Estimate species richness in a community----
out2 <- raw %>% count(Scientific.Name) %>% 
  column_to_rownames(var = "Scientific.Name") %>% 
  as.data.frame() %>%
  SpadeR::ChaoSpecies(datatype = "abundance", k=3, conf = 0.95) #k is cut-ff point for rare species
out2$Basic_data_information
out2$Species_table

#Environmental variable on species----
## Regression model----
### GLM for continous variable----

#Calculate the abundance of the species Treron oxyurus by transect
species_abundance <- raw %>%
  filter(Scientific.Name == "Treron oxyurus") %>% #change the species interest
  group_by(Transect) %>%
  summarise(Abundance = n())

#join dataset with environment variable
raw_env_bird <- left_join(raw_env, species_abundance, by = "Transect") %>%
  replace_na(list(Abundance = 0))

glm1 <- glm(Abundance ~ Elev, data = raw_env_bird, family = "poisson") #fit a GLM with Poisson family
glm2 <- glm(Abundance ~ Elev, data = raw_env_bird, family = "quasipoisson") #fit a GLM with Quasi-Poisson family

#check both result
summary(glm1)
summary(glm2)

#for quasi-poisson, the model is better when dispersion is not close to 1, so we choose glm1

# Predict abundance based on the best model
pred_df <- data.frame(Elev = seq(min(raw_env_bird$Elev), max(raw_env_bird$Elev), length.out = 100))
pred_df$Abundance <- predict(glm1, newdata = pred_df, type = "response")

# Create a ggplot with the raw_env_bird data and the predicted Abundance values
ggplot(raw_env_bird, aes(x = Elev, y = Abundance)) +
  geom_point(alpha = 0.5, color = "#3B6B9A") + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, color = "#8C510A") +
  theme_minimal(base_size = 16) +
  ggtitle("Abundance vs. Elevation") +
  xlab("Elevation") +
  ylab("Abundance") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = "bottom")

### GLM for categorical variable----
glm3 <- glm(Abundance ~ Tuplah, data = raw_env_bird, family = "poisson")
glm4 <- glm(Abundance ~ Tuplah, data = raw_env_bird, family = "quasipoisson")

#for quasi-poisson, the model is better when dispersion is not close to 1, so we choose glm4

# Generate predicted values and standard errors for each level of Tuplah
pred_df <- data.frame(Tuplah = unique(raw_env_bird$Tuplah))
pred_df$Predicted <- predict(glm4, newdata = pred_df, type = "response", se.fit = TRUE)$fit
pred_df$SE <- predict(glm4, newdata = pred_df, type = "response", se.fit = TRUE)$se.fit

# Create a ggplot with the raw_env_bird data and the predicted Abundance values
ggplot(data = pred_df, aes(x = Tuplah, y = Predicted)) +
  geom_bar(stat = "identity", fill = "#6BAED6") +
  geom_errorbar(aes(ymin = Predicted - SE, ymax = Predicted + SE), width = 0.2, color = "#3182BD") +
  labs(x = "Landcover", y = "Predicted count",
       title = "Predicted counts of abundance by landcover",
       fill = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


##CCA (Canonical Correspondence Analysis)----

# Create abundance matrix with species data
sp.abd <- raw %>%
  count(Transect, Scientific.Name) %>%
  pivot_wider(names_from = Scientific.Name, values_from = n, values_fill = list(n = 0)) %>%
  column_to_rownames(var = "Transect")

# Create abundance matrix with species data
sp.env <- raw_env %>%
  column_to_rownames(var = "Transect") %>%
  select (c(3:7)) %>%
  mutate(Elev = scale(Elev)) %>%
  mutate(Tree.count = scale(Tree.count)) %>%
  mutate(Imp.Tree = scale(Imp.Tree)) %>%
  mutate(Tuplah = as.factor(Tuplah))

#Calculate the Canonical Correspondence Analysis
OM1 <- cca(sp.abd ~ Elev + Tuplah + Tree.count, sp.env)

# Generate a simple plot with biplots, species scores and site scores
plot1 <- plot(OM1, display=c("bp", "sp", "sites"))

# generate fancier plot

# Convert the CCA model to a data frame
df <- fortify(OM1)

#make plot
ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), name = "CCA1") +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), name = "CCA2") +    
  geom_point(aes(x=CCA1, y=CCA2), data=filter(df, Score=="sites"), size = 3, color = "#1f77b4") + 
  geom_segment(aes(x = 0, y = 0, xend = CCA1*2.5, yend = CCA2*2.5), data = filter(df, Score=="biplot"),color = "#d62728", arrow=arrow(length = unit(0.3, "cm"))) +
  geom_text_repel(aes(x=CCA1*2.5, y=CCA2*2.5, label=Label),data = filter(df, Score=="biplot"),color = "#d62728", size = 4) + 
  geom_text_repel(aes(x=CCA1, y=CCA2, label=Label), data=filter(df, Score=="species"),color = "#2ca02c", size = 4) +
  geom_text_repel(aes(x=CCA1, y=CCA2, label=Label), data=filter(df, Score=="sites"),color = "#1f77b4", size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")