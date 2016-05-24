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
##  Univariate normal datasets  ##
##################################

# Generate univariate normal datasets.

N <- array(data = list(NULL), dim = 15)

N[[1]]$density <- "Gaussian"
N[[1]]$w <- 1
N[[1]]$Theta <- list(pdf1 = "normal", theta1.1 = 0, theta2.1 = 1)

N[[2]]$density <- "Skewed unimodal"
N[[2]]$w <- c(rep(1/5, 2), 3/5)
N[[2]]$Theta <- list(pdf1 = "normal", theta1.1 = c(0, 1/2, 13/12), theta2.1 = c(1, 2/3, 5/9))

N[[3]]$density <- "Strongly skewed"
N[[3]]$w <- rep(1/8, 8)

theta1 <- function(x) 3 * ((2/3)^x - 1)
theta2 <- function(x) (2/3)^x

N[[3]]$Theta <- list(pdf1 = "normal", theta1.1 = theta1(0:7), theta2.1 = theta2(0:7))

N[[4]]$density <- "Kurtotic unimodal"
N[[4]]$w <- c(2/3, 1/3)
N[[4]]$Theta <- list(pdf1 = "normal", theta1.1 = rep(0, 2), theta2.1 = c(1, 1/10))

N[[5]]$density <- "Outlier"
N[[5]]$w <- c(1/10, 9/10)
N[[5]]$Theta <- list(pdf1 = "normal", theta1.1 = rep(0, 2), theta2.1 = c(1, 1/10))

N[[6]]$density <- "Bimodal"
N[[6]]$w <- rep(1/2, 2)
N[[6]]$Theta <- list(pdf1 = "normal", theta1.1 = c(-1, 1), theta2.1 = rep(2/3, 2))

N[[7]]$density <- "Separated bimodal"
N[[7]]$w <- rep(1/2, 2)
N[[7]]$Theta <- list(pdf1 = "normal", theta1.1 = c(-3/2, 3/2), theta2.1 = rep(1/2, 2))

N[[8]]$density <- "Skewed bimodal"
N[[8]]$w <- c(3/4, 1/4)
N[[8]]$Theta <- list(pdf1 = "normal", theta1.1 = c(0, 3/2), theta2.1 = c(1, 1/3))

N[[9]]$density <- "Trimodal"
N[[9]]$w <- c(rep(9/20, 2), 1/10)
N[[9]]$Theta <- list(pdf1 = "normal", theta1.1 = c(-6/5, 6/5, 0), theta2.1 = c(rep(3/5, 2), 1/4))

N[[10]]$density <- "Claw"
N[[10]]$w <- c(1/2, rep(1/10, 5))

theta1 <- function(x) x / 2 - 1

N[[10]]$Theta <- list(pdf1 = "normal", theta1.1 = c(0, theta1(0:4)), theta2.1 = c(1, rep(1/10, 5)))

N[[11]]$density <- "Double claw"
N[[11]]$w <- c(rep(49/100, 2), rep(1/350, 7))

theta1 <- function(x) (x - 3) / 2

N[[11]]$Theta <- list(pdf1 = "normal", theta1.1 = c(-1, 1, theta1(0:6)), theta2.1 = c(rep(2/3, 2), rep(1/100, 7)))

N[[12]]$density <- "Asimetric claw"

w <- function(x) 2^(1 - x) / 31

N[[12]]$w <- c(1/2, w(-2:2))

theta1 <- function(x) x + 1 / 2
theta2 <- function(x) 2^(-x) / 10

N[[12]]$Theta <- list(pdf1 = "normal", theta1.1 = c(0, theta1(-2:2)), theta2 = c(1, theta2(-2:2)))

N[[13]]$density <- "Asimetric double claw"
N[[13]]$w <- c(rep(46/100, 2), rep(1/300, 3), rep(7/300, 3))

theta1.1 <- function(x) 2 * x - 1 
theta1.2 <- function(x) -x / 2 
theta1.3 <- function(x) x / 2 

N[[13]]$Theta <- list(pdf1 = "normal", theta1.1 = c(theta1.1(0:1), theta1.2(1:3), theta1.3(1:3)), theta2.1 = c(rep(2/3, 2), rep(1/100, 3), rep(7/100, 3)))

N[[14]]$density <- "Smooth comb"

w <- function(x) 2^(5 - x) / 63

N[[14]]$w <- w(0:5)

theta1 <- function(x) (65 - 96 * (1 / 2)^x) / 21
theta2 <- function(x) (32 / 63) / 2^x

N[[14]]$Theta <- list(pdf1 = "normal", theta1.1 = theta1(0:5), theta2.1 = theta2(0:5))

N[[15]]$density <- "Discrete comb"
N[[15]]$w <- c(rep(2/7, 3), rep(1/21, 3))

theta1.1 <- function(x) (12 * x - 15) / 7
theta1.2 <- function(x) 2 * x / 7

N[[15]]$Theta <- list(pdf1 = "normal", theta1.1 = c(theta1.1(0:2), theta1.2(8:10)), theta2.1 = c(rep(2/7, 3), rep(1/21, 3)))

# Set dataset size.

n <- c(100, 1000, 10000)

normal <- array(data = list(NULL), dim = c(15, 3))

for (i in 1:15) {
  for (j in 1:length(n)) {
    normal[[i, j]] <- RNGMIX(Dataset.name = paste("normal_", i, "_", j, "_", 1:100, sep = ""), 
      n = ceiling(N[[i]]$w * n[j]), Theta = N[[i]]$Theta)
  }
}

# Estimate number of components, component weights and component parameters.

normalest <- array(data = list(NULL), dim = c(15, 3))
table <- array(data = 0, dim = c(15, 6))

for (j in 1:length(n)) {
  Sturges <- as.integer(1 + log2(n[j])) # Minimum v follows the Sturges rule.
  Log10 <- as.integer(10 * log10(n[j])) # Maximum v follows the Log10 rule.
  RootN <- as.integer(2 * n[j]^0.5) # Maximum v follows the RootN rule.

  for (i in 1:15) {
    normalest[[i, j]] <- REBMIX(Dataset = normal[[i, j]]@Dataset,
      Preprocessing = "histogram", cmax = 20, Criterion = "BIC",
      pdf = "normal", K = Sturges:Log10, y0 = 0.0, 
      ymin = normal[[i, j]]@ymin, ymax = normal[[i, j]]@ymax)

    table[i, j * 2 - 1] <- mean(as.numeric(normalest[[i, j]]@summary$c))  
    table[i, j * 2] <- sd(as.numeric(normalest[[i, j]]@summary$c))
  }
}

table
