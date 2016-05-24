### R code from vignette source 'mRMRe.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
options(keep.source=TRUE)


###################################################
### code chunk number 2: install-pkg (eval = FALSE)
###################################################
## install.packages("mRMRe")


###################################################
### code chunk number 3: loadlib
###################################################
library(mRMRe)


###################################################
### code chunk number 4: utils
###################################################
set.thread.count(2)


###################################################
### code chunk number 5: loadlib
###################################################
data(cgps)
data.annot <- data.frame(cgps.annot)
data.cgps <- data.frame(cgps.ic50, cgps.ge)


###################################################
### code chunk number 6: mim
###################################################
## Test on a dummy dataset.

# Create a dummy data set
library(survival)
df <- data.frame(
    "surv1" = Surv(runif(100),
                   sample(0:1, 100, replace = TRUE)),
    "cont1" = runif(100),
    "disc1" = factor(sample(1:5, 100, replace = TRUE),
                     ordered = TRUE),
    "surv2" = Surv(runif(100),
                   sample(0:1, 100, replace = TRUE)),
    "cont2" = runif(100),
    "cont3" = runif(100),
    "surv3" = Surv(runif(100),
                   sample(0:1, 100, replace = TRUE)),
    "disc2" = factor(sample(1:5, 100, replace = TRUE),
                     ordered = TRUE))
dd <- mRMR.data(data = df)

# Show a partial mutual information matrix.
print(mim(subsetData(dd, 1:4, 1:4)))


###################################################
### code chunk number 7: mim2
###################################################
## Test on the 'cgps' dataset, where the
## variables are all of continuous type.

dd <- mRMR.data(data = data.cgps)
dd <- subsetData(dd, 1:10, 1:10)

# Uses Spearman as correlation estimator
spearman_mim <- mim(dd, continuous_estimator = "spearman") 
print(spearman_mim[1:4, 1:4])

# Uses Pearson as correlation estimator
pearson_mim <- mim(dd, continuous_estimator = "pearson") 
print(pearson_mim[1:4, 1:4])


###################################################
### code chunk number 8: correlations
###################################################
# Compute c-index between feature 1 and 2
correlate(cgps.ge[, 1], cgps.ge[, 2], method = "cindex")

# Compute Cramer's V
x <- sample(factor(c("CAT_1", "CAT_2", "CAT_3"),
                   ordered = TRUE), 100, replace = TRUE)
y <- sample(factor(c("CAT_1", "CAT_2"),
                   ordered = TRUE), 100, replace = TRUE)
correlate(x, y, method = "cramersv")

# Compute Pearson coefficient with random strata and
# sample weights between features 1 and 2
strata <- sample(factor(c("STRATUM_1", "STRATUM_2",
                          "STRATUM_3"),
                       ordered = TRUE), 
	               nrow(cgps.ge), replace = TRUE)
weights <- runif(nrow(cgps.ge))
correlate(cgps.ge[, 1], cgps.ge[, 2], strata = strata,
          weights = weights, method = "pearson")


###################################################
### code chunk number 9: classic.mRMR
###################################################
dd <- mRMR.data(data = data.cgps)

mRMR.classic(data = dd, target_indices = c(1),
             feature_count = 30)


###################################################
### code chunk number 10: ensemble.mRMR
###################################################
dd <- mRMR.data(data = data.cgps)

# For mRMR.classic-like results
mRMR.ensemble(data = dd, target_indices = c(1),
              solution_count = 1, feature_count = 30)

# For mRMR.ensemble-like results
mRMR.ensemble(data = dd, target_indices = c(1),
              solution_count = 5, feature_count = 30)


###################################################
### code chunk number 11: causality
###################################################
ensemble <- mRMR.ensemble(data = dd, target_indices = c(1),
                          solution_count = 5,
                          feature_count = 10)
causality(ensemble)


###################################################
### code chunk number 12: sessionInfo
###################################################
toLatex(sessionInfo())


