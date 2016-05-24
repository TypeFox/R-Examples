##############################################
## R sources for reproducing the results in ##
##   Marko Nagode:                          ##
##   Finite Mixture Modeling via REBMIX     ##
##############################################

options(prompt = "> ", continue = "+ ", width = 70,
  useFancyQuotes = FALSE, digits = 3)

###################
## Preliminaries ##
###################

# load package.

library("rebmix")

##################################
## Multivariate normal datasets ##
##################################

# Generate multivariate normal datasets.

n <- c(75, 100, 125, 150, 175)

Theta <- list(pdf1 = rep("normal", 4),
  theta1.1 = c(10, 12, 10, 12),
  theta2.1 = c(1, 1, 1, 1),
  pdf2 = rep("normal", 4),
  theta1.2 = c(8.5, 10.5, 8.5, 10.5),
  theta2.2 = c(1, 1, 1, 1),
  pdf3 = rep("normal", 4),
  theta1.3 = c(12, 14, 12, 14),
  theta2.3 = c(1, 1, 1, 1),
  pdf4 = rep("normal", 4),
  theta1.4 = c(13, 15, 7, 9),
  theta2.4 = c(2, 2, 2, 2),
  pdf5 = rep("normal", 4),
  theta1.5 = c(7, 9, 13, 15),
  theta2.5 = c(3, 3, 3, 3))

normal <- RNGMIX(Dataset = paste("normal_", 1:100, sep = ""), n = n, Theta = Theta)

# Estimate number of components, component weights and component parameters.

Sturges <- as.integer(1 + log2(sum(n))) # Minimum v follows the Sturges rule.
Log10 <- as.integer(10 * log10(sum(n))) # Maximum v follows the Log10 rule.

normalest <- REBMIX(Dataset = normal@Dataset, Preprocessing = "histogram",
  K = Sturges:Log10, Criterion = "BIC", pdf = rep("normal", 4))

c <- as.numeric(normalest@summary$c)
IC <- as.numeric(normalest@summary$IC)

summary(c)
summary(IC, digits = 5)

format(length(c[c == 5]) / length(c))
