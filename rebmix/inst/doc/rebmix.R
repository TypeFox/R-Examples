### R code from vignette source 'rebmix.Rnw'
### Encoding: CP1250

###################################################
### code chunk number 1: rebmix-code
###################################################
##############################################
## R sources for reproducing the results in ##
##   Marko Nagode:                          ##
##   rebmix: The Rebmix Package             ##
##############################################

options(prompt = "R> ", continue = "+  ", width = 80,
  useFancyQuotes = FALSE, digits = 3)


###################################################
### code chunk number 2: rebmix-code
###################################################
###################
## Preliminaries ##
###################

## load package and set prompt before starting new page to TRUE.

library("rebmix")
devAskNewPage(ask = TRUE)


###################################################
### code chunk number 3: rebmix-code
###################################################
######################
##  Gamma datasets  ##
######################

## Generate gamma datasets.

n <- c(100, 100, 100, 100)

Theta <- list(pdf1 = "gamma",
  theta1.1 = c(1/100, 1/100, 1/100, 1/100),
  theta2.1 = c(200, 400, 600, 800))

gamma1 <- RNGMIX(Dataset.name = "gamma1", n = n, Theta = Theta)

n <- c(40, 360)

Theta <- list(pdf1 = "gamma",
  theta1.1 = c(1/27, 1/270),
  theta2.1 = c(9, 90))

gamma2 <- RNGMIX(Dataset.name = "gamma2", n = n, Theta = Theta)

n <- c(80, 240, 80)

Theta <- list(pdf1 = "gamma",
  theta1.1 = c(1/20, 1, 1/20),
  theta2.1 = c(40, 6, 200))

gamma3 <- RNGMIX(Dataset.name = "gamma3", n = n, Theta = Theta)


###################################################
### code chunk number 4: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

gamma1est <- REBMIX(Dataset = gamma1@Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = c("AIC", "BIC"),
  pdf = "gamma",
  K = 30:80)

gamma2est <- REBMIX(Dataset = gamma2@Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  K = 30:80)

gamma3est <- REBMIX(Dataset = gamma3@Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  K = 30:80)


###################################################
### code chunk number 5: rebmix-code
###################################################
summary(gamma2est)

coef(gamma3est)


###################################################
### code chunk number 6: rebmix-code
###################################################
## Bootstrap finite mixture.
gamma3boot <- boot(x = gamma3est, pos = 1, Bootstrap = "p", B = 10)

gamma3boot

summary(gamma3boot)


###################################################
### code chunk number 7: gamma2-fig
###################################################
plot(gamma2est, pos = 1, what = c("den", "dis"), ncol = 2, npts = 1000)


###################################################
### code chunk number 8: rebmix-code
###################################################
#########################
##   Poisson dataset   ##
#########################

## Generate the Poisson dataset.

n <- c(200, 200, 200)

Theta <- list(pdf1 = rep("Poisson", 2),
  theta1.1 = c(3, 2),
  theta2.1 = c(NA, NA),
  pdf2 = rep("Poisson", 2),
  theta1.2 = c(9, 10),
  theta2.2 = c(NA, NA),
  pdf3 = rep("Poisson", 2),
  theta1.3 = c(15, 16),
  theta2.3 = c(NA, NA))

poisson <- RNGMIX(Dataset.name = paste("Poisson_", 1:10, sep = ""), n = n, Theta = Theta)


###################################################
### code chunk number 9: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

poissonest <- REBMIX(Dataset = poisson@Dataset,
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "MDL5",
  pdf = rep("Poisson", 2),
  K = 1)


###################################################
### code chunk number 10: rebmix-code
###################################################
## Visualize results.

summary(poissonest)

coef(poissonest, pos = 9)


###################################################
### code chunk number 11: poisson-fig
###################################################
plot(poissonest, pos = 7, what = c("dens", "marg", "IC", "D", "logL"), nrow = 2, ncol = 3, npts = 1000)


###################################################
### code chunk number 12: poisson-clu-fig
###################################################
poissonclu <- RCLRMIX(x = poissonest, pos = 9, Zt = poisson@Zt)

plot(poissonclu)


###################################################
### code chunk number 13: rebmix-code
###################################################
data("wreath", package = "mclust")

## Estimate number of components, component weights and component parameters.

n <- nrow(wreath)

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows the Sturges rule.
  as.integer(2 * sum(n)^0.5)) # Maximum v follows the RootN rule.

wreathest <- REBMIX(model = "REBMVNORM",
  Dataset = list(as.data.frame(wreath)),
  Preprocessing = "histogram",
  cmax = 20,
  Criterion = "BIC",
  pdf = rep("normal", ncol(wreath)),
  K = K[1]:K[2])


###################################################
### code chunk number 14: rebmix-code
###################################################
summary(wreathest)

coef(wreathest)


###################################################
### code chunk number 15: wreath-fig
###################################################
plot(wreathest)


###################################################
### code chunk number 16: wreath-clu-fig
###################################################
wreathclu <- RCLRMIX(model = "RCLRMVNORM", x = wreathest)

plot(wreathclu)


###################################################
### code chunk number 17: rebmix-code
###################################################
data("Baudry_etal_2010_JCGS_examples", package = "mclust")

## Estimate number of components, component weights and component parameters.

n <- nrow(ex4.1)

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows the Sturges rule.
  as.integer(2 * sum(n)^0.5)) # Maximum v follows the RootN rule.

ex4.1est <- REBMIX(model = "REBMVNORM",
  Dataset = list(as.data.frame(ex4.1)),
  Preprocessing = "Parzen window",
  cmax = 10,
  Criterion = "AIC",
  pdf = rep("normal", ncol(ex4.1)),
  K = K[1]:K[2])


###################################################
### code chunk number 18: rebmix-code
###################################################
summary(ex4.1est)


###################################################
### code chunk number 19: ex4_1-fig
###################################################
plot(ex4.1est, pos = 1, what = c("dens"), nrow = 1, ncol = 1)


###################################################
### code chunk number 20: ex4_1-clu-fig
###################################################
ex4.1clu <- RCLRMIX(model = "RCLRMVNORM", x = ex4.1est)

plot(ex4.1clu)


###################################################
### code chunk number 21: rebmix-code
###################################################
data("iris")

# Show level attributes discrete variables.

levels(iris[["Class"]])

# Split iris dataset into three subsets for three Classes
# and remove Class column.

iris_set <- subset(iris, subset = Class == "iris-setosa", select = c(-Class))
iris_ver <- subset(iris, subset = Class == "iris-versicolor", select = c(-Class))
iris_vir <- subset(iris, subset = Class == "iris-virginica", select = c(-Class))


###################################################
### code chunk number 22: rebmix-code
###################################################
# Split datasets into train (75%) and test (25%) subsets.

set.seed(5)

Prob <- 0.75

n_set <- nrow(iris_set); s_set <- sample.int(n = n_set, size = as.integer(n_set * Prob))

iris_set_train <- iris_set[s_set,]; iris_set_test <- iris_set[-s_set,]

n_ver <- nrow(iris_ver); s_ver <- sample.int(n = n_ver, size = as.integer(n_ver * Prob))

iris_ver_train <- iris_ver[s_ver,]; iris_ver_test <- iris_ver[-s_ver,]

n_vir <- nrow(iris_vir); s_vir <- sample.int(n = n_vir, size = as.integer(n_vir * Prob))

iris_vir_train <- iris_vir[s_vir,]; iris_vir_test <- iris_vir[-s_vir,]

iris_test = rbind(iris_set_test, iris_ver_test, iris_vir_test)


###################################################
### code chunk number 23: rebmix-code
###################################################
Zt <- factor(c(rep(0, nrow(iris_set_test)),
  rep(1, nrow(iris_ver_test)),
  rep(2, nrow(iris_vir_test))))


###################################################
### code chunk number 24: rebmix-code
###################################################
# Estimate number of components, component weights and component
# parameters for train subsets.

n <- range(nrow(iris_set_train), nrow(iris_ver_train), nrow(iris_vir_train))

K <- c(as.integer(1 + log2(sum(n[1]))), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n[2]))) # Maximum v follows log10 rule.

K <- c(floor(K[1]^(1/4)), ceiling(K[2]^(1/4)))

irisest <- REBMIX(model = "REBMVNORM",
  Dataset = list(iris_set_train = iris_set_train,
                 iris_ver_train = iris_ver_train,
                 iris_vir_train = iris_vir_train),
  Preprocessing = "Parzen window",
  cmax = 10,
  Criterion = "ICL-BIC",
  pdf = rep("normal", 4),
  K = K[1]:K[2])


###################################################
### code chunk number 25: rebmix-code
###################################################
iriscla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(irisest),
  Dataset = iris_test,
  Zt = Zt)


###################################################
### code chunk number 26: rebmix-code
###################################################
iriscla

summary(iriscla)


###################################################
### code chunk number 27: iris-cla-fig
###################################################
plot(iriscla, nrow = 3, ncol = 2)


###################################################
### code chunk number 28: rebmix-code
###################################################
data("adult")

# Find complete cases.

adult <- adult[complete.cases(adult),]

# Replace levels with numbers.

adult <- as.data.frame(data.matrix(adult))


###################################################
### code chunk number 29: rebmix-code
###################################################
# Split adult dataset into two train subsets for two Incomes
# and remove Type and Income columns.

trainle50k <- subset(adult, subset = (Type == 2) & (Income == 1),
  select = c(-Type, -Income))
traingt50k <- subset(adult, subset = (Type == 2) & (Income == 2),
  select = c(-Type, -Income))

trainall <- subset(adult, subset = Type == 2, select = c(-Type, -Income))

train <- as.factor(subset(adult, subset = Type == 2, select = c(Income))[, 1])


###################################################
### code chunk number 30: rebmix-code
###################################################
# Extract test dataset form adult dataset and remove Type
# and Income columns.

testle50k <- subset(adult, subset = (Type == 1) & (Income == 1),
  select = c(-Type, -Income))
testgt50k <- subset(adult, subset = (Type == 1) & (Income == 2),
  select = c(-Type, -Income))

testall <- subset(adult, subset = Type == 1, select = c(-Type, -Income))

test <- as.factor(subset(adult, subset = Type == 1, select = c(Income))[, 1])


###################################################
### code chunk number 31: rebmix-code
###################################################
# Estimate number of components, component weights and component
# parameters for Naive Bayes.

cmax <- unlist(lapply(apply(trainall, 2, unique), length))

adultest <- list(0)

for (i in 1:14) {
  adultest[[i]] <- REBMIX(Dataset = list(as.data.frame(trainle50k[, i]),
    as.data.frame(traingt50k[, i])),
    Preprocessing = "histogram",
    cmax = if (cmax[i] > 120) 12 else cmax[i],
    Criterion = "BIC",
    pdf = if (cmax[i] > 120) "normal" else "Dirac",
    K = if (cmax[i] > 120) 13:43 else 1)
}


###################################################
### code chunk number 32: rebmix-code
###################################################
# Best-first feature subset selection.

c <- NULL; rvs <- 1:14; Error <- 1.0

for (i in 1:14) {
  k <- NA

  for (j in rvs) {
    adultcla <- RCLSMIX(x = adultest[c(c, j)],
      Dataset = as.data.frame(trainall[, c(c, j)]),
      Zt = train)

    if (adultcla@Error < Error) {
      Error <- adultcla@Error; k <- j
    }
  }

  if (is.na(k)) {
    break
  }
  else {
    c <- c(c, k); rvs <- rvs[-which(rvs == k)]
  }
}

# Error on train dataset.

Error


###################################################
### code chunk number 33: rebmix-code
###################################################
# Selected features.

adultcla <- RCLSMIX(x = adultest[c],
  Dataset = as.data.frame(testall[, c]),
  Zt = test)


###################################################
### code chunk number 34: rebmix-code
###################################################
adultcla

summary(adultcla)


###################################################
### code chunk number 35: adult-cla-fig
###################################################
plot(adultcla, nrow = 5, ncol = 2)


###################################################
### code chunk number 36: rebmix-code
###################################################
rm(list = ls())


