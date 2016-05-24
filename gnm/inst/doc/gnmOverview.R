### R code from vignette source 'gnmOverview.Rnw'

###################################################
### code chunk number 1: gnmOverview.Rnw:65-66
###################################################
getOption("SweaveHooks")[["eval"]]()
options(SweaveHooks = list(eval = function() options(show.signif.stars = FALSE)))


###################################################
### code chunk number 2: Load_gnm
###################################################
getOption("SweaveHooks")[["eval"]]()
library(gnm)


###################################################
### code chunk number 3: migrationData
###################################################
getOption("SweaveHooks")[["eval"]]()
count <- c(11607,   100,   366,   124,
              87, 13677,   515,   302,
             172,   225, 17819,   270,
              63,   176,   286, 10192 )
region <- c("NE", "MW", "S", "W")
row <-  gl(4, 4, labels = region)
col <-  gl(4, 1, length = 16, labels = region)


###################################################
### code chunk number 4: squareTableModels
###################################################
getOption("SweaveHooks")[["eval"]]()
independence <- glm(count ~ row + col, family = poisson)
quasi.indep <- glm(count ~ row + col + Diag(row, col), family = poisson)
symmetry <- glm(count ~ Symm(row, col), family = poisson)
quasi.symm <- glm(count ~ row + col + Symm(row, col), family = poisson)
comparison1 <- anova(independence, quasi.indep, quasi.symm)
print(comparison1, digits = 7)
comparison2 <- anova(symmetry, quasi.symm)
print(comparison2)


###################################################
### code chunk number 5: EriksonData
###################################################
getOption("SweaveHooks")[["eval"]]()
### Collapse to 7 by 7 table as in Erikson et al. (1982)
erikson <- as.data.frame(erikson)
lvl <- levels(erikson$origin)
levels(erikson$origin) <- levels(erikson$destination) <-
    c(rep(paste(lvl[1:2], collapse = " + "), 2), lvl[3],
      rep(paste(lvl[4:5], collapse = " + "), 2), lvl[6:9])
erikson <- xtabs(Freq ~ origin + destination + country, data = erikson)


###################################################
### code chunk number 6: wedderburn
###################################################
getOption("SweaveHooks")[["eval"]]()
##  data from Wedderburn (1974), see ?barley
logitModel <- glm(y ~ site + variety, family = wedderburn, data = barley)
fit <- fitted(logitModel)
print(sum((barley$y - fit)^2 / (fit * (1-fit))^2))


###################################################
### code chunk number 7: termPredictors
###################################################
getOption("SweaveHooks")[["eval"]]()
print(temp <- termPredictors(quasi.symm))
rowSums(temp) - quasi.symm$linear.predictors


###################################################
### code chunk number 8: RC_homogeneous_model_1
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
RChomog1 <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus, verbose = FALSE)


###################################################
### code chunk number 9: RC_homogeneous_model_2
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(2)
RChomog2 <- update(RChomog1)


###################################################
### code chunk number 10: Compare_coefficients
###################################################
getOption("SweaveHooks")[["eval"]]()
compareCoef <- cbind(coef(RChomog1), coef(RChomog2))
colnames(compareCoef) <- c("RChomog1", "RChomog2")
round(compareCoef, 4)


###################################################
### code chunk number 11: Summarize_model
###################################################
getOption("SweaveHooks")[["eval"]]()
summary(RChomog2)


###################################################
### code chunk number 12: RC_homogeneous_constrained_model1
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
RChomogConstrained1 <- update(RChomog1, constrain = length(coef(RChomog1)))


###################################################
### code chunk number 13: RC_homogeneous_constrained_model2
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(2)
RChomogConstrained2 <- update(RChomogConstrained1)
identical(coef(RChomogConstrained1), coef(RChomogConstrained2))


###################################################
### code chunk number 14: Eliminate_Eg
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
n <- 1000
x <- rep(rnorm(n), rep(3, n))
counts <- as.vector(rmultinom(n, 10, c(0.7, 0.1, 0.2)))
rowID <- gl(n, 3, 3 * n)
resp <- gl(3, 1, 3 * n)


###################################################
### code chunk number 15: Double_UNIDIFF_model
###################################################
getOption("SweaveHooks")[["eval"]]()
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
                     + Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote), family = poisson,
                     data = cautres)


###################################################
### code chunk number 16: Contrast_matrix
###################################################
getOption("SweaveHooks")[["eval"]]()
coefs <- names(coef(doubleUnidiff))
contrCoefs <- coefs[grep(", religion:vote", coefs)]
nContr <- length(contrCoefs)
contrMatrix <- matrix(0, length(coefs), nContr,
                      dimnames = list(coefs, contrCoefs))
contr <- contr.sum(contrCoefs)
# switch round to contrast with first level
contr <- rbind(contr[nContr, ], contr[-nContr, ])
contrMatrix[contrCoefs, 2:nContr] <- contr
contrMatrix[contrCoefs, 2:nContr]


###################################################
### code chunk number 17: Check_estimability_1
###################################################
getOption("SweaveHooks")[["eval"]]()
checkEstimable(doubleUnidiff, contrMatrix)


###################################################
### code chunk number 18: Check_estimability_2
###################################################
getOption("SweaveHooks")[["eval"]]()
coefs <- names(coef(doubleUnidiff))
contrCoefs <- coefs[grep("[.]religion", coefs)]
nContr <- length(contrCoefs)
contrMatrix <- matrix(0, length(coefs), length(contrCoefs),
                      dimnames = list(coefs, contrCoefs))
contr <- contr.sum(contrCoefs)
contrMatrix[contrCoefs, 2:nContr] <- rbind(contr[nContr, ], contr[-nContr, ])
checkEstimable(doubleUnidiff, contrMatrix)


###################################################
### code chunk number 19: Get_contrasts_1
###################################################
getOption("SweaveHooks")[["eval"]]()
myContrasts <- getContrasts(doubleUnidiff,
                            pickCoef(doubleUnidiff, ", religion:vote"))
myContrasts


###################################################
### code chunk number 20: qvplot
###################################################
getOption("SweaveHooks")[["eval"]]()
plot(myContrasts,
  main = "Relative strength of religion-vote association, log scale",
xlab = "Election", levelNames = 1:4)


###################################################
### code chunk number 21: RCmodel
###################################################
getOption("SweaveHooks")[["eval"]]()
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS + Mult(SES, MHS),
                family = poisson, data = mentalHealth)


###################################################
### code chunk number 22: RCmodel_constrained
###################################################
getOption("SweaveHooks")[["eval"]]()
RC1model2 <- gnm(count ~ SES + MHS + Mult(1, SES, MHS),
                 constrain = "[.]SES[AF]", constrainTo = c(0, 1),
                 ofInterest = "[.]SES",
                 family = poisson, data = mentalHealth)
summary(RC1model2)


###################################################
### code chunk number 23: getContrasts_simple
###################################################
getOption("SweaveHooks")[["eval"]]()
getContrasts(RC1model, pickCoef(RC1model, "[.]SES"), ref = "first",
             scaleRef = "first", scaleWeights = c(rep(0, 5), 1))


###################################################
### code chunk number 24: two-way
###################################################
getOption("SweaveHooks")[["eval"]]()
xtabs(y ~ site + variety, barley)


###################################################
### code chunk number 25: residSVD
###################################################
getOption("SweaveHooks")[["eval"]]()
emptyModel <- gnm(y ~ -1, family = wedderburn, data = barley)
biplotStart <- residSVD(emptyModel, barley$site, barley$variety, d = 2)
biplotModel <- gnm(y ~ -1 + instances(Mult(site, variety), 2),
family = wedderburn, data = barley, start =  biplotStart)


###################################################
### code chunk number 26: residSVDplot
###################################################
getOption("SweaveHooks")[["eval"]]()
plot(coef(biplotModel), biplotStart,
     main = "Comparison of residSVD and MLE for a 2-dimensional
 biplot model", ylim = c(-2, 2), xlim = c(-4, 4))
abline(a = 0, b = 1, lty = 2)


###################################################
### code chunk number 27: Set_contrasts_attribute
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)


###################################################
### code chunk number 28: RC1_model
###################################################
getOption("SweaveHooks")[["eval"]]()
RC1model <- gnm(count ~ SES + MHS + Mult(SES, MHS), family = poisson,
                data = mentalHealth)
RC1model


###################################################
### code chunk number 29: Normalize_scores
###################################################
getOption("SweaveHooks")[["eval"]]()
rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
rowScores <- coef(RC1model)[10:15]
colScores <- coef(RC1model)[16:19]
rowScores <- rowScores - sum(rowScores * rowProbs)
colScores <- colScores - sum(colScores * colProbs)
beta1 <- sqrt(sum(rowScores^2 * rowProbs))
beta2 <- sqrt(sum(colScores^2 * colProbs))
assoc <- list(beta = beta1 * beta2,
              mu = rowScores / beta1,
              nu = colScores / beta2)
assoc


###################################################
### code chunk number 30: Elliptical_contrasts
###################################################
getOption("SweaveHooks")[["eval"]]()
mu <- getContrasts(RC1model, pickCoef(RC1model, "[.]SES"),
                   ref = rowProbs, scaleWeights = rowProbs)
nu <- getContrasts(RC1model, pickCoef(RC1model, "[.]MHS"),
                   ref = colProbs, scaleWeights = colProbs)
mu
nu


###################################################
### code chunk number 31: RC2_model
###################################################
getOption("SweaveHooks")[["eval"]]()
RC2model <- gnm(count ~ SES + MHS + instances(Mult(SES, MHS), 2),
                family = poisson, data = mentalHealth)
RC2model


###################################################
### code chunk number 32: Homogeneous_effects
###################################################
getOption("SweaveHooks")[["eval"]]()
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)
RChomog


###################################################
### code chunk number 33: Heterogeneous_effects
###################################################
getOption("SweaveHooks")[["eval"]]()
RCheterog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Mult(origin, destination), family = poisson,
               data = occupationalStatus)
anova(RChomog, RCheterog)


###################################################
### code chunk number 34: Transform_to_counts
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)


###################################################
### code chunk number 35: Class_mobility
###################################################
getOption("SweaveHooks")[["eval"]]()
classMobility <- gnm(yvar ~ Dref(origin, destination),
                       family = binomial, data = voting)
classMobility


###################################################
### code chunk number 36: Class_mobility_weights
###################################################
getOption("SweaveHooks")[["eval"]]()
DrefWeights(classMobility)


###################################################
### code chunk number 37: Salariat_factors
###################################################
getOption("SweaveHooks")[["eval"]]()
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)


###################################################
### code chunk number 38: Social_mobility
###################################################
getOption("SweaveHooks")[["eval"]]()
socialMobility <- gnm(yvar ~ Dref(origin, destination,
                                  delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)
socialMobility


###################################################
### code chunk number 39: social_mobility_weights
###################################################
getOption("SweaveHooks")[["eval"]]()
DrefWeights(socialMobility)


###################################################
### code chunk number 40: Downward_mobility
###################################################
getOption("SweaveHooks")[["eval"]]()
downwardMobility <- gnm(yvar ~ Dref(origin, destination,
                                    delta = ~ 1 + downward),
                        family = binomial, data = voting)
downwardMobility
DrefWeights(downwardMobility)


###################################################
### code chunk number 41: UNIDIFF_model
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
unidiff <- gnm(Freq ~ educ*orig + educ*dest + Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", family = poisson,
               data = yaish, subset = (dest != 7))
coef(unidiff)


###################################################
### code chunk number 42: Unidiff_contrasts
###################################################
getOption("SweaveHooks")[["eval"]]()
getContrasts(unidiff, ofInterest(unidiff))


###################################################
### code chunk number 43: double_UNIDIFF_model
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                     Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote),

             family = poisson, data = cautres)
getContrasts(doubleUnidiff, rev(pickCoef(doubleUnidiff, ", class:vote")))
getContrasts(doubleUnidiff, rev(pickCoef(doubleUnidiff, ", religion:vote")))


###################################################
### code chunk number 44: Scale_yields
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- interaction(wheat$tillage, wheat$summerCrop, wheat$manure,
                         wheat$N, sep = "")


###################################################
### code chunk number 45: AMMI_model
###################################################
getOption("SweaveHooks")[["eval"]]()
mainEffects <- gnm(yield.scaled ~ year + treatment, family = gaussian,
                   data = wheat)
svdStart <- residSVD(mainEffects, year, treatment, 3)
bilinear1 <- update(mainEffects, . ~ . + Mult(year, treatment),
                    start = c(coef(mainEffects), svdStart[,1]))


###################################################
### code chunk number 46: AOD
###################################################
getOption("SweaveHooks")[["eval"]]()
anova(mainEffects, bilinear1, test = "F")


###################################################
### code chunk number 47: AMMI_model2
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
barleyModel <- gnm(height ~ year + genotype + Mult(year, genotype),
                   data = barleyHeights)


###################################################
### code chunk number 48: Spherical_contrasts
###################################################
getOption("SweaveHooks")[["eval"]]()
gamma <- getContrasts(barleyModel, pickCoef(barleyModel, "[.]y"),
                      ref = "mean", scaleWeights = "unit")
delta <- getContrasts(barleyModel, pickCoef(barleyModel, "[.]g"),
                      ref = "mean", scaleWeights = "unit")
gamma
delta


###################################################
### code chunk number 49: CI
###################################################
getOption("SweaveHooks")[["eval"]]()
gamma[[2]][,1] + (gamma[[2]][,2]) %o% c(-1.96, 1.96)
delta[[2]][,1] + (delta[[2]][,2]) %o% c(-1.96, 1.96)


###################################################
### code chunk number 50: SVD
###################################################
getOption("SweaveHooks")[["eval"]]()
svd(termPredictors(barleyModel)[, "Mult(year, genotype)"])$d


###################################################
### code chunk number 51: Biplot_model
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(83)
biplotModel <- gnm(y ~ -1 + instances(Mult(site, variety), 2),
                   family = wedderburn, data = barley)


###################################################
### code chunk number 52: Row_and_column_scores
###################################################
getOption("SweaveHooks")[["eval"]]()
barleyMatrix <- xtabs(biplotModel$predictors ~ site + variety,
                      data = barley)
barleySVD <- svd(barleyMatrix)
A <- sweep(barleySVD$u, 2, sqrt(barleySVD$d), "*")[, 1:2]
B <- sweep(barleySVD$v, 2, sqrt(barleySVD$d), "*")[, 1:2]
rownames(A) <- levels(barley$site)
rownames(B) <- levels(barley$variety)
colnames(A) <- colnames(B) <- paste("Component", 1:2)
A
B


###################################################
### code chunk number 53: Biplot1
###################################################
getOption("SweaveHooks")[["eval"]]()
barleyCol <- c("red", "blue")
plot(rbind(A, B), pch = c(levels(barley$site), levels(barley$variety)),
     col = rep(barleyCol, c(nlevels(barley$site), nlevels(barley$variety))),
     xlim = c(-4, 4), ylim = c(-4, 4), main = "Biplot for barley data",
     xlab = "Component 1", ylab = "Component 2")
text(c(-3.5, -3.5), c(3.9, 3.6), c("sites: A-I","varieties: 1-9, X"),
     col = barleyCol, adj = 0)


###################################################
### code chunk number 54: Biplot2
###################################################
getOption("SweaveHooks")[["eval"]]()
plot(rbind(A, B), pch = c(levels(barley$site), levels(barley$variety)),
     col = rep(barleyCol, c(nlevels(barley$site), nlevels(barley$variety))),
     xlim = c(-4, 4), ylim = c(-4, 4), main = "Biplot for barley data",
     xlab = "Component 1", ylab = "Component 2")
text(c(-3.5, -3.5), c(3.9, 3.6), c("sites: A-I","varieties: 1-9, X"),
     col = barleyCol, adj = 0)
abline(a = 0, b = tan(pi/3))
abline(a = 0, b = -tan(pi/6))
abline(a = 2.6, b = tan(pi/3), lty = 2)
abline(a = 4.5, b = tan(pi/3), lty = 2)
abline(a = 1.3, b = -tan(pi/6), lty = 2)
text(2.8, 3.9, "v-axis", font = 3)
text(3.8, -2.7, "h-axis", font = 3)


###################################################
### code chunk number 55: Double_additive
###################################################
getOption("SweaveHooks")[["eval"]]()
variety.binary <- factor(match(barley$variety, c(2,3,6), nomatch = 0) > 0,
                        labels = c("rest", "2,3,6"))
doubleAdditive <- gnm(y ~ variety + Mult(site, variety.binary),
                      family = wedderburn, data = barley)


###################################################
### code chunk number 56: Compare_chi-squared
###################################################
getOption("SweaveHooks")[["eval"]]()
biplotModChiSq <- sum(residuals(biplotModel, type = "pearson")^2)
doubleAddChiSq <- sum(residuals(doubleAdditive, type = "pearson")^2)
c(doubleAddChiSq - biplotModChiSq,
  doubleAdditive$df.residual - biplotModel$df.residual)


###################################################
### code chunk number 57: Re-express_data
###################################################
getOption("SweaveHooks")[["eval"]]()
set.seed(1)
subset(backPain, x1 == 1 & x2 == 1 & x3 == 1)
backPainLong <- expandCategorical(backPain, "pain")
head(backPainLong)


###################################################
### code chunk number 58: Stereotype_model
###################################################
getOption("SweaveHooks")[["eval"]]()
oneDimensional <- gnm(count ~ pain + Mult(pain, x1 + x2 + x3),
                      eliminate = id, family = "poisson", data = backPainLong)
oneDimensional


###################################################
### code chunk number 59: Qualitative_model
###################################################
getOption("SweaveHooks")[["eval"]]()
threeDimensional  <- gnm(count ~ pain + pain:(x1 + x2 + x3), eliminate = id,
                         family = "poisson", data = backPainLong)
threeDimensional


###################################################
### code chunk number 60: Calculate_log-likelihood
###################################################
getOption("SweaveHooks")[["eval"]]()
logLikMultinom <- function(model, size){
    object <- get(model)
    l <- sum(object$y * log(object$fitted/size))
    c(nParameters = object$rank - nlevels(object$eliminate), logLikelihood = l)
}
size <- tapply(backPainLong$count, backPainLong$id, sum)[backPainLong$id]
t(sapply(c("oneDimensional", "threeDimensional"), logLikMultinom, size))


###################################################
### code chunk number 61: Constrain_slopes
###################################################
getOption("SweaveHooks")[["eval"]]()
## before constraint
summary(oneDimensional)
oneDimensional <- gnm(count ~ pain + Mult(pain, offset(x1) + x2 + x3),
                      eliminate = id, family = "poisson", data = backPainLong)
## after constraint
summary(oneDimensional)


###################################################
### code chunk number 62: Get_slopes
###################################################
getOption("SweaveHooks")[["eval"]]()
getContrasts(oneDimensional, pickCoef(oneDimensional, "[.]pain"))


###################################################
### code chunk number 63: singleExp
###################################################
getOption("SweaveHooks")[["eval"]]()
x <- 1:100
y <- exp(- x / 10)
set.seed(1)
saved.fits <- list()
for (i in 1:100) saved.fits[[i]] <- gnm(y ~ Exp(1 + x), verbose = FALSE)
table(zapsmall(sapply(saved.fits, deviance)))


###################################################
### code chunk number 64: singleExp2
###################################################
getOption("SweaveHooks")[["eval"]]()
saved.fits[[2]]


###################################################
### code chunk number 65: doubleExp
###################################################
getOption("SweaveHooks")[["eval"]]()
x <- 1:100
y <- exp(- x / 10) + 2 * exp(- x / 50)
set.seed(1)
saved.fits <- list()
for (i in 1:100) {
    saved.fits[[i]] <- suppressWarnings(gnm(y ~ Exp(1 + x, inst = 1) +
                                            Exp(1 + x, inst = 2),
                                            verbose = FALSE))
}
table(round(unlist(sapply(saved.fits, deviance)), 4))


###################################################
### code chunk number 66: doubleExp2
###################################################
getOption("SweaveHooks")[["eval"]]()
singleExp <- gnm(y ~ Exp(1 + x), start = c(NA, NA, -0.1), verbose = FALSE)
singleExp
meanOnly <- gnm(y ~ 1, verbose = FALSE)
meanOnly
plot(x, y, main = "Two sub-optimal fits to a sum-of-exponentials curve")
lines(x, fitted(singleExp))
lines(x, fitted(meanOnly), lty = "dashed")


###################################################
### code chunk number 67: doubleExp3
###################################################
getOption("SweaveHooks")[["eval"]]()
gnm(y ~ instances(Exp(1 + x), 2), start = c(NA, NA, -0.1, NA, -0.1),
    verbose = FALSE)


