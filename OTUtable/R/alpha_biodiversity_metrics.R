# Calculations for alpha biodiversity

# Observed richness = number of taxa observed in a sample
obs_richness <- function(sample){
  richness <- length(which(sample > 0))
  return(richness)
}

# Shannon Index: H = -sum(Pi * ln(Pi)) where Pi is each species in a sample in proportional form
shannon <- function(sample){
  sample <- sample[which(sample > 0)]
  prop <- sample / sum(sample)
  H <- -sum(prop * log(prop))
  return(H)
}

# Pielou's Evenness: E = Shannon/ln(# taxa observed)
pielou <- function(sample){
  sample <- sample[which(sample > 0)]
  prop <- sample / sum(sample)
  H <- -sum(prop * log(prop))
  Hmax <- log(length(sample))
  return(H / Hmax)
}

# Chao1 Richness = observed richness + (# of singletons^2)/(2 * # of doubletons)
chao1 <- function(sample){
  sample <- sample[which(sample > 0)]
  obs <- length(sample)
  single <- length(which(sample == 1))
  double <- length(which(sample == 2))
  richness <- obs + (single ^ 2)/(2 * double)
  return(richness) 
}

