## ------------------------------------------------------------------------
# Set a seed for reproducibility
set.seed(123) 

# Define subject types.

# Truthfully respond "Yes" to direct question
N.trueadmitter <- 500

# Falsely respond "No" to direct question
N.withholder <- 500

# Truthfully respond "No" to direct question
N.innocent <- 500

type <- rep(c("TA", "WH", "IN"), times=c(N.trueadmitter, N.withholder, N.innocent))

## ------------------------------------------------------------------------
D <- ifelse(type=="TA", 1, 0)
direct.est <- mean(D)
direct.est

## ------------------------------------------------------------------------
N <- length(type)
# Generate list response potential outcomes

# Control potential outcome
Y0 <- sample(1:4, N, replace=TRUE)

# Treated potential outcome is 1 higher for true admitters and withholders
Y1 <- Y0 + ifelse(type %in% c("TA", "WH"), 1, 0)

# Conduct random assignment
Z <- rbinom(N, 1, 0.5)

# Reveal list responses
Y <- Z*Y1 + (1-Z)*Y0

list.est <- mean(Y[Z==1]) - mean(Y[Z==0])
list.se <- sqrt((var(Y[Z==1])/sum(Z) + var(Y[Z==0])/sum(1-Z)))
list.est
list.se

## ------------------------------------------------------------------------
library(list)
# Wrap up all data in a dataframe
df <- data.frame(Y, Z, D)
out.1 <- combinedListDirect(formula = Y ~ Z, data = df, treat = "Z", direct = "D")
out.1

## ------------------------------------------------------------------------
summary(out.1)

## ------------------------------------------------------------------------
# Define three subject types as before plus one new type

N.trueadmitter <- 400
N.withholder <- 500
N.innocent <- 500

# Truthfully responds "Yes" to direct question
# but decreases response to the non-sensitive items 
# in the presence of the sensitive item
N.designaffected <- 100

type <- rep(c("TA", "WH", "IN", "DA"), 
            times=c(N.trueadmitter, N.withholder, N.innocent, N.designaffected))
N <- length(type)

D <- ifelse(type%in%c("TA","DA"), 1, 0)

# Control potential outcome
Y0 <- sample(1:4, N, replace=TRUE)

# Treated potential outcome is 1 higher for true admitters and withholders
# Note that it is NOT higher for those who are "design affected"
Y1 <- Y0 + ifelse(type %in% c("TA", "WH"), 1, 0)

Z <- rbinom(N, 1, 0.5)
Y <- Z*Y1 + (1-Z)*Y0
df <- data.frame(Y, Z, D)

out.2 <- combinedListDirect(formula = Y ~ Z, data = df, treat = "Z", direct = "D")

# Extract Placebo Test I results 
unlist(out.2$placebo.I)

## ------------------------------------------------------------------------
# Define three subject types as before plus one new type

N.trueadmitter <- 400
N.withholder <- 500
N.innocent <- 500

# Truthfully answers "Yes" when in control
# But falsely answers "No" when in treatment
N.affectedbytreatment <- 100

type <- rep(c("TA", "WH", "IN", "ABT"), 
            times=c(N.trueadmitter, N.withholder, N.innocent, N.affectedbytreatment))
N <- length(type)

# Direct Question Potential outcomes
D0 <- ifelse(type%in%c("TA","ABT"), 1, 0)
D1 <- ifelse(type%in%c("TA"), 1, 0)

# List Experiment potential outcomes
Y0 <- sample(1:4, N, replace=TRUE)
Y1 <- Y0 + ifelse(type %in% c("TA", "WH"), 1, 0)

# Reveal outcomes according to random assignment
Z <- rbinom(N, 1, 0.5)
Y <- Z*Y1 + (1-Z)*Y0
D <- Z*D1 + (1-Z)*D0

df <- data.frame(Y, Z, D)

out.3 <- combinedListDirect(formula = Y ~ Z, data = df, treat = "Z", direct = "D")

# Extract Placebo Test II results 
unlist(out.3$placebo.II)

## ------------------------------------------------------------------------
# Define subject types.
N.trueadmitter <- 500
N.withholder <- 500
N.innocent <- 500

type <- rep(c("TA", "WH", "IN"), times=c(N.trueadmitter, N.withholder, N.innocent))
N <- length(type)

# Generate a predictive pre-treatment covariate "X")
X <- rnorm(N, sd = 2)

# Control potential outcome is related to "X"
Y0 <- as.numeric(cut(X + runif(N), breaks = 4))
Y1 <- Y0 + ifelse(type %in% c("TA", "WH"), 1, 0)

Z <- rbinom(N, 1, 0.5)
D <- ifelse(type=="TA", 1, 0)
Y <- Z*Y1 + (1-Z)*Y0

df <- data.frame(Y, Z, D, X)

# Conduct estimation without covariate adjustment
out.4 <- combinedListDirect(formula = Y ~ Z, data = df, treat = "Z", direct = "D")
out.4

# Conduct estimation with covariate adjustment
# Just add the covariate on the right-hand side of the formula
out.5 <- combinedListDirect(formula = Y ~ Z + X, data = df, treat = "Z", direct = "D")
out.5

