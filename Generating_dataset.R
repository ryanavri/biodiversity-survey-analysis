##Generate data set for birds survey----
#THIS CHAPTER IS OPTIONAL!!, YOU MAY SKIP TO THE ANLYSIS OF DIVERSITY.
#it just gave you an example how to produce "mock" dataset from GBIF.
#the prepared dataset from "data.RData" contain;
#raw for bird dataset
#raw_m for mammal dataset
#raw_t for trees dataset
#raw_env for environmental dataset 

#Set probabilities for reproducible code
set.seed(123) #pseudo random number
n <- 666 #number of sample

# Set probability of selected value to make the observation for each species and transect uneven
probs_a <- c(rep(0.1, 2), rep(0.05, 4), rep(0.005, 3), rep(0.0025, 1))
probs_b <- runif(90, min = 0.001, max = 0.03)
probs_b[c(1:5, 8:10)] <- 0.9 # make some species more dominant
probs_b[c(6, 11)] <- 0.0001 # make some species less selected

# Pick selected taxon from GBIF dataset
#birds
b <- GGPNP_GBIF %>%
  filter(class=="Aves") %>% 
  select(species) %>%
  unique()
b <- b[sample(nrow(b), 90), ]

#trees
t <- GGPNP_GBIF %>%
  filter(class=="Magnoliopsida") %>% 
  select(species) %>%
  unique()
t <- t[sample(nrow(b), 90), ]

# Create dataframe for bird
raw <- data.frame(
  Transect = sample(c("TS01", "TS02", "TS03", "TS04", "TS05", "TS06","TS07", "TS08", "TS09", "TS10"), n, replace = TRUE, prob = probs_a),
  Scientific.Name = sample(unique(b$Scientific.Name), n, replace = TRUE, prob = probs_b),
  Distance = round(runif(n, min = 0, max = 30),0))

# Create dataframe for tree
raw.t <- data.frame(
  Transect = sample(c("TS01", "TS02", "TS03", "TS04", "TS05", "TS06","TS07", "TS08", "TS09", "TS10"), n, replace = TRUE, prob = probs_a),
  Scientific.Name = sample(unique(t$species), n, replace = TRUE, prob = probs_b),
  dbh = round(runif(n, min = 5, max = 80),0))


# Distance sampling----
## Generate dataset for primates----
# Set the number of observations
n <- 80

# Set the maximum distance
max_dist <- 30

# Generate a vector of distances
distances <- seq(0, max_dist, length.out = n)

# Generate detection probabilities based on the assumption
detection_probs <- exp(-0.1 * distances) + 0.5 * exp(-0.2 * (distances - 15)^2)

# Define the species names
species_names <- c("Macaca fascicularis", "Trachypithecus auratus", "Presbytis comata", "Hylobates javanicus", "Nycticebus javanicus")

# Define the detection probability multipliers for the last two species
multipliers <- c(1, 1, 0.5, 0.3, 0.1)  # Adjust the values as desired

# Simulate observed distances with detection probabilities for each transect and species
simulated_data <- data.frame()
for (i in 1:10) {
  transect <- paste0("TS", sprintf("%02d", i))
  species <- sample(species_names, size = n, replace = TRUE, prob = multipliers)
  observed_distances <- sample(distances, size = n, replace = TRUE, prob = detection_probs)
  transect_data <- data.frame(Transect = transect, Species = species, Distance = observed_distances)
  simulated_data <- rbind(simulated_data, transect_data)
}

simdat_primate_ds <- simulated_data