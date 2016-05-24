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

# Load package.

library("rebmix")

#######################################
## Mixed continuous-discrete dataset ##
#######################################

# Generate mixed continuous-discrete dataset.

n <- c(400, 100, 500)

Theta <- list(pdf1 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.1 = c(1, 2, 10, 2), 
  theta2.1 = c(0.3, NA, 0.9, 3),
  pdf2 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.2 = c(3.5, 10, 10, 10),
  theta2.2 = c(0.2, NA, 0.1, 7),
  pdf3 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.3 = c(2.5, 15, 10, 25), 
  theta2.3 = c(0.4, NA, 0.7, 20))

mixed <- RNGMIX(Dataset = "mixed", n = n, Theta = Theta)

# Estimate number of components, component weights and component parameters.

Sturges <- as.integer(1 + log2(sum(n))) # Minimum v follows the Sturges rule.
Log10 <- as.integer(10 * log10(sum(n))) # Maximum v follows the Log10 rule.

mixedest <- REBMIX(Dataset = mixed@Dataset, Preprocessing = "histogram", cmax = 9,
  Criterion = "BIC", pdf = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1 = c(NA, NA, 10, NA), K = kseq(Sturges, Log10, 0.1))

#library("tikzDevice") # Uncomment to use tikzDevice package.
#tikz("mixed.tex", width = 4.5, height = 4.5) # Uncomment to use tikzDevice package.
plot(mixedest, what = c("dens", "marg", "IC", "logL"), nrow = 4, ncol = 3, npts = 200)
#dev.off() # Uncomment to use tikzDevice package.

# Bootstrap finite mixture.

mixedboot <- boot(x = mixedest, pos = 1, Bootstrap = "p", B = 100)

mixedboot
