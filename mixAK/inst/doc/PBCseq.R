###
###  Full pdf document describing the code included here is available at
###  http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf
###
### ==============================================================================


###################################################
### code chunk number 1: Options
###################################################
OPT <- options(width = 135, digits = 3)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.1                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 4: Load package and data
###################################################
library("mixAK")
data("PBC910", package = "mixAK")
tail(PBC910)[, c(1, 2, 3, 6:9)]


###################################################
### code chunk number 5: Show additional info for data
###################################################
dim(PBC910)
head(PBC910)
summary(PBC910)
length(unique(PBC910[, "id"]))


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.2                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 6: Extract longitudinal profiles and plot profiles of the variable lbili
###################################################
ip <- getProfiles(t = "month", 
                  y = c("lbili", "platelet", "spiders", "jspiders"), 
                  id = "id", data = PBC910)
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
             main = "Log(bilirubin)", highlight = c(1, length(ip)),
             xlab = "Time (months)", ylab = "Log(bilirubin)")


###################################################
### code chunk number 7: Print extracted longitudinal profiles of the first and last subject
###################################################
print(ip[[1]])
print(ip[[length(ip)]])


###################################################
### code chunk number 9: 01-long_prof
###################################################
par(mar = c(4, 4, 4, 0) + 0.1, bty = "n")
iShow <- c(1, length(ip))     ### indeces of subjects to be highlighted

layout(autolayout(3))
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
       auto.layout = FALSE, main = "Log(bilirubin)", 
       xlab = "Time (months)", ylab = "Log(bilirubin)",
       highlight = iShow)
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
       auto.layout = FALSE, main = "Log(platelet)", 
       trans = log, xlab = "Time (months)", ylab = "Log(platelet)",
       highlight = iShow)
plotProfiles(ip = ip, data = PBC910, var = "jspiders", tvar = "month", 
       lines = FALSE, points = TRUE,
       auto.layout = FALSE, main = "Blood vessel malform.", 
       xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)",
       highlight = iShow)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.4                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 10: Running MCMC (eval = FALSE)
###################################################
set.seed(20042007)
mod <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = PBC910[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = PBC910[, "month"]),
    z = list(lbili    = PBC910[, "month"], 
             platelet = PBC910[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    prior.b = list(Kmax = 2), 
    nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),
    parallel = FALSE)


###################################################
### code chunk number 11: Class and elements of object mod
###################################################
class(mod)
names(mod)


###################################################
### code chunk number 12: Class and elements of the first element of mod
###################################################
class(mod[[1]])
names(mod[[1]])


###################################################
### code chunk number 13: order and rank objects
### (discussed only in an extended version of the paper)
###################################################
mod[[1]]$order_b[1:3,]
mod[[1]]$rank_b[1:3,]


###################################################
### code chunk number 14: Running relabeling algorithm and storing of sampled component probabilities (eval = FALSE)
###################################################
mod <- NMixRelabel(mod, type = "stephens", keep.comp.prob = TRUE)


###################################################
### code chunk number 15: Parallel running of relabeling algorithm (eval = FALSE)
###################################################
## mod <- NMixRelabel(mod, type = "stephens", keep.comp.prob = TRUE, parallel = TRUE)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.5                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 17: First ten sampled values of the model deviance
###################################################
print(mod[[1]]$Deviance[1:10], digits = 9)


###################################################
### code chunk number 18: First ten sampled values of the fixed effects
###################################################
print(mod[[1]]$alpha[1:10,])


###################################################
### code chunk number 19: First ten sampled values of the square roots of the dispersion parameter
###################################################
print(mod[[1]]$sigma_eps[1:10,])


###################################################
### code chunk number 20: First ten sampled values of the overall mean and variance parameters of the random effects
###################################################
print(mod[[1]]$mixture_b[1:10,])


###################################################
### code chunk number 21: Autocorrelation for deviance chains
###################################################
library("coda")
DevChains <- mcmc.list(mcmc(mod[[1]]$Deviance), mcmc(mod[[2]]$Deviance))
autocorr(DevChains)


###################################################
### code chunk number 22: Relabeled sample of the mixture means
###################################################
muSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "mu_b")
print(muSamp1[1:3,])


###################################################
### code chunk number 23: Relabeled sample of the mixture weights
###################################################
wSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "w_b")
print(wSamp1[1:3,])


###################################################
### code chunk number 24: Relabeled sample of the mixture variances
###################################################
varSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "var_b")
print(varSamp1[1:3,])


###################################################
### code chunk number 25: Relabeled sample of the mixture standard deviations
###################################################
sdSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "sd_b")
print(sdSamp1[1:3,])


###################################################
### code chunk number 26: Relabeled sample of the mixture correlations
###################################################
corSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "cor_b")
print(corSamp1[1:3,])


###################################################
### code chunk number 27: Relabeled sample of the mixture covariance matrices - lower triangles
###################################################
DSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "Sigma_b")
print(DSamp1[1:3,])


###################################################
### code chunk number 28: Relabeled sample of the mixture inverted covariance matrices - lower triangles
###################################################
QSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "Q_b")
print(QSamp1[1:3,])


###################################################
### code chunk number 29: Relabeled sample of the Cholesky factors of the mixture inverted covariance matrices - lower triangles
###################################################
LiSamp1 <- NMixChainComp(mod[[1]], relabel = TRUE, param = "Li_b")
print(LiSamp1[1:3,])


###################################################
### code chunk number 31: 02-trace_deviance
###################################################
par(mar = c(4, 4, 0, 1) + 0.1)
tracePlots(mod, param = "Deviance")


###################################################
### code chunk number 32: Traceplot of fixed effects and residual standard deviation
###################################################
tracePlots(mod, param = "alpha")
tracePlots(mod, param = "sigma_eps")


###################################################
### code chunk number 33: Traceplot of Eb and standard deviations and correlations derived from varb
###################################################
tracePlots(mod, param = "Eb")
tracePlots(mod, param = "SDb")
tracePlots(mod, param = "Corb")


###################################################
### code chunk number 34: traceplots for mixture weighs and means and standard deviations
###################################################
tracePlots(mod, param = "w_b")
tracePlots(mod, param = "mu_b")
tracePlots(mod, param = "sd_b")


###################################################
### code chunk number 35: Traceplots for mixture weighs and means and standard deviations after relabeling
###################################################
tracePlots(mod, param = "w_b",  relabel = TRUE)
tracePlots(mod, param = "mu_b", relabel = TRUE)
tracePlots(mod, param = "sd_b", relabel = TRUE)


###################################################
### code chunk number 36: traceplots of variance hyperparameters (eval = FALSE)
###################################################
tracePlots(mod, param = "gammaInv_b")
tracePlots(mod, param = "gammaInv_eps")


###################################################
### code chunk number 38: 03-trace_Emub_chain2
###################################################
par(mar = c(4, 4, 0, 1) + 0.1)
tracePlots(mod, param = "mu_b", relabel = TRUE)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.6                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 39: Posterior summary statistics for basic model parameters
###################################################
print(mod)


###################################################
### code chunk number 40: Coda posterior summary for regression parameters
###################################################
name.Eb <- paste("b.Mean.", 1:5, sep = "")
RegrChain1 <- cbind(mod[[1]]$mixture_b[, name.Eb], 
                    mod[[1]]$alpha, mod[[1]]$sigma_eps)
RegrChain2 <- cbind(mod[[2]]$mixture_b[, name.Eb], 
                    mod[[2]]$alpha, mod[[2]]$sigma_eps)
colnames(RegrChain1) <- colnames(RegrChain2) <- 
    c(paste(rep(c("lbili", "platelet", "spiders"), each = 2), 
            ":", rep(c("Intcpt", "Slope"), 3), sep=""),
      "lbili:res_std_dev")
RegrChain1 <- mcmc(RegrChain1)
RegrChain2 <- mcmc(RegrChain2)
summary(mcmc.list(RegrChain1, RegrChain2))


###################################################
### code chunk number 41: Coda posterior summary for SDb and Corb
###################################################
name.SDb <- paste("b.SD.", 1:5, sep = "")
SDbChain1 <- mcmc(mod[[1]]$mixture_b[, name.SDb]) 
SDbChain2 <- mcmc(mod[[2]]$mixture_b[, name.SDb]) 
summary(mcmc.list(SDbChain1, SDbChain2))
#
name.Corb <- paste("b.Corr.", c(2:5, 3:5, 4:5, 5), ".", rep(1:4, 4:1), sep = "")
CorbChain1 <- mcmc(mod[[1]]$mixture_b[, name.Corb])
CorbChain2 <- mcmc(mod[[2]]$mixture_b[, name.Corb])
summary(mcmc.list(CorbChain1, CorbChain2))


###################################################
### code chunk number 42: HPD intervals for regression parameters
###################################################
HPDinterval(mcmc.list(RegrChain1, RegrChain2))


###################################################
### code chunk number 43: HPD intervals for SDb and Corb
###################################################
HPDinterval(mcmc.list(SDbChain1, SDbChain2))
HPDinterval(mcmc.list(CorbChain1, CorbChain2))


###################################################
### code chunk number 44: Posterior densities of beta (eval = FALSE)
###################################################
COL <- rep(rainbow_hcl(3, start = 30, end = 210), each = 2)
par(mfcol = c(2, 3))
for (i in 1:6){
  densplot(RegrChain1[, i], show.obs = FALSE, col = COL[i], lwd = 2)
  title(main = colnames(RegrChain1)[i])
}


###################################################
### code chunk number 45: Posterior means of mixture parameters
###################################################
NMixSummComp(mod[[1]])


###################################################
### code chunk number 46: Coda summary for the mixture means based on the first sampled chain
###################################################
summary(mcmc(muSamp1))
HPDinterval(mcmc(muSamp1))


###################################################
### code chunk number 47: Coda summary for the mixture weights based on the first sampled chain
###################################################
summary(mcmc(wSamp1))
HPDinterval(mcmc(wSamp1))


###################################################
### code chunk number 48: Coda summary for the mixture variance based on the first sampled chain
###################################################
summary(mcmc(varSamp1))
HPDinterval(mcmc(varSamp1))


###################################################
### code chunk number 49: Coda summary for the mixture standard deviations based on the first sampled chain
###################################################
summary(mcmc(sdSamp1))
HPDinterval(mcmc(sdSamp1))


###################################################
### code chunk number 50: Coda summary for the mixture correlations based on the first sampled chain
###################################################
summary(mcmc(corSamp1))
HPDinterval(mcmc(corSamp1))


###################################################
### code chunk number 51: Coda summary for the mixture covariance matrices based on the first sampled chain
###################################################
summary(mcmc(DSamp1))
HPDinterval(mcmc(DSamp1))


###################################################
### code chunk number 52: Coda summary for the mixture inverted covariance matrices based on the first sampled chain
###################################################
summary(mcmc(QSamp1))
HPDinterval(mcmc(QSamp1))


###################################################
### code chunk number 53: Coda summary for the Cholesky decompositions of the mixture inverted covariance matrices based on the first sampled chain
###################################################
summary(mcmc(LiSamp1))
HPDinterval(mcmc(LiSamp1))


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.7                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 54: Calculation of cluster specific fitted longitudinal profiles based on posterior means using Gaussian quadrature (eval = FALSE)
###################################################
delta <- 0.3
tpred <- seq(0, 30, by = delta)
fit <- fitted(mod[[1]], x =list("empty", "empty", tpred), 
                        z =list(tpred, tpred, "empty"),
              glmer = TRUE)
names(fit) <- c("lbili", "platelet", "spiders")


###################################################
### code chunk number 56: Inaccurate calculation of cluster specific fitted longitudinal profiles based on posterior means without Gaussian quadrature
###################################################
fit0 <- fitted(mod[[1]], x = list("empty", "empty", tpred), 
                         z = list(tpred, tpred, "empty"),
               glmer = FALSE)
names(fit0) <- c("lbili", "platelet", "spiders")


###################################################
### code chunk number 57: Print part of calculated longitudinal profiles for platelet counts
###################################################
print(fit[["platelet"]][1:3,], digits = 5)


###################################################
### code chunk number 59: 04-fitted_prof
###################################################
par(mar = c(4, 4, 4, 0) + 0.1, bty = "n")
K <- mod[[1]]$prior.b$Kmax
clCOL <- c("darkgreen", "red3")
obsCOL <- "azure3"

layout(autolayout(3))
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
       auto.layout = FALSE, main = "Log(bilirubin)", 
       xlab = "Time (months)", ylab = "Log(bilirubin)", col = obsCOL)
for (k in 1:K) lines(tpred, fit[["lbili"]][, k], col = clCOL[k], lwd = 2)
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
       auto.layout = FALSE, main = "Platelet count", 
       xlab = "Time (months)", ylab = "Platelet count", col = obsCOL)
for (k in 1:K) lines(tpred, fit[["platelet"]][, k], col = clCOL[k], lwd = 2)
plotProfiles(ip = ip, data = PBC910, var = "jspiders", tvar = "month", 
       lines = FALSE, points = TRUE,
       auto.layout = FALSE, main = "Blood vessel malform.", 
       xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)", 
       bg = obsCOL, col = "grey70")
for (k in 1:K) lines(tpred, fit[["spiders"]][, k], col = clCOL[k], lwd = 2)


###################################################
### code chunk number 60: Posterior sample of the component probabilities
###################################################
print(mod[[1]]$comp.prob[1:3, 1:8])


###################################################
### code chunk number 61: Print posterior means of the ICPs for the first three subjects
###################################################
print(mod[[1]]$poster.comp.prob[1:3,])


###################################################
### code chunk number 62: Classification based on posterior means of the individual component probabilities
###################################################
groupMean <- apply(mod[[1]]$poster.comp.prob, 1, which.max)
pMean <- apply(mod[[1]]$poster.comp.prob, 1, max)
table(groupMean)


###################################################
### code chunk number 63: Print posterior medians of the individual component probabilities for the first three subjects
###################################################
print(mod[[1]]$quant.comp.prob[["50%"]][1:3,])


###################################################
### code chunk number 64: Print posterior 2.5 and 97.5 percent quantiles of the individual component probabilities for the first three subjects
###################################################
print(mod[[1]]$quant.comp.prob[["2.5%"]][1:3,])
print(mod[[1]]$quant.comp.prob[["97.5%"]][1:3,])


###################################################
### code chunk number 65: Classification based on posterior medians of the individual component probabilities
###################################################
groupMed <- apply(mod[[1]]$quant.comp.prob[["50%"]], 1, which.max)
pMed <- apply(mod[[1]]$quant.comp.prob[["50%"]], 1, max)
table(groupMed)
table(groupMean, groupMed)


###################################################
### code chunk number 66: Subjects with different classification using posterior mean and median
###################################################
pMeanMed <- data.frame(Mean = pMean, Median = pMed)
rownames(pMeanMed) <- unique(PBC910$id)
print(pMeanMed[groupMean != groupMed, ])


###################################################
### code chunk number 68: 05-post_dist_p
###################################################
IDS <- unique(PBC910$id)

N <- ncol(mod[[1]]$comp.prob) / K
ID <- c(2, 7, 11)

par(mfrow = c(1, 3), mar = c(4, 4, 4, 1) + 0.1, bty = "n")
for (id in ID){
  i <- (1:N)[IDS == id]
  hist(mod[[1]]$comp.prob[, (i - 1) * K + 1], xlim = c(0, 1), prob = TRUE, 
     xlab = expression(paste("P(U=1|Y=y; ", theta, ")", 
                             sep = "")), 
     col = rainbow_hcl(1, start=60), 
     main = paste("ID ", id, "   (", 
       format(round(pMean[i], 3), nsmall = 3), ",  ", 
       format(round(pMed[i], 3), nsmall = 3), ")", 
       sep=""))  
}  


###################################################
### code chunk number 69: Classification using the HPD credible intervals of the individual component probabilities
###################################################
pHPD <- HPDinterval(mcmc(mod[[1]]$comp.prob))
pHPDlower <- matrix(pHPD[, "lower"], ncol = 2, byrow = TRUE)
pHPDupper <- matrix(pHPD[, "upper"], ncol = 2, byrow = TRUE)
rownames(pHPDlower) <- rownames(pHPDupper) <- unique(PBC910$id)
groupHPD <- groupMed
groupHPD[groupHPD == 1 & pHPDlower[, 1] <= 0.5] <- NA
groupHPD[groupHPD == 2 & pHPDlower[, 2] <= 0.5] <- NA
table(groupHPD, useNA = "ifany")


###################################################
### code chunk number 70: Variable groupMeanHPD
### (this is discussed only in an extended version of the paper)
###################################################
groupMeanHPD <- as.character(groupMean)
groupMeanHPD[is.na(groupHPD) & groupMean == 1] <- "1_NA"
groupMeanHPD[is.na(groupHPD) & groupMean == 2] <- "2_NA"
groupMeanHPD <- factor(groupMeanHPD, levels = c("1", "2", "1_NA", "2_NA"))
table(groupMeanHPD)


###################################################
### code chunk number 71: Add group indicators to data
###################################################
TAB <- table(PBC910$id)
PBC910$groupMed <- factor(rep(groupMed, TAB))
PBC910$groupHPD <- factor(rep(groupHPD, TAB))
PBC910$groupHPD <- addNA(PBC910$groupHPD)


###################################################
### code chunk number 72: Extract longitudinal profiles including the group indicators
###################################################
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "jspiders", "groupMed", "groupHPD"),
   id = "id", data = PBC910)
print(ip[[1]])


###################################################
### code chunk number 75: 06-long_prof_by_group
###################################################
oldPar <- par(mar = c(5, 4, 3, 0) + 0.1, bty = "n", cex.main = 1.1, cex.lab = 0.9, cex.axis = 0.9)
GCOL <- rainbow_hcl(3, start = 220, end = 40, c = 50, l = 60)[c(2, 3, 1)]
GCOL2 <- c("darkgreen", "red3", "skyblue")
names(GCOL) <- names(GCOL2) <- levels(PBC910$groupHPD)

layout(autolayout(4))

# Log(bilirubin):
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
   gvar = "groupHPD", col = GCOL, 
   auto.layout = FALSE, main = "Log(bilirubin)",
   xlab = "Time (months)", ylab = "Log(bilirubin)")

# Legend:
plot(c(0, 100), c(0, 100), type = "n", 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(0, 90, legend = c("Group 1", "Group 2", 
                         "Uncertain classification"),
       lty = 1, lwd = 5, col = GCOL, y.intersp = 1.5, cex=1.1)

# Platelet count:
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
   gvar = "groupHPD", col = GCOL, 
   auto.layout = FALSE, main = "Platelet count",
   xlab = "Time (months)", ylab = "Platelet count")

# Blood vessel malformations
plotProfiles(ip = ip, data = PBC910, var = "jspiders",  tvar = "month",
   lines = FALSE, points = TRUE,
   gvar = "groupHPD", bg = GCOL, col = GCOL2, 
   auto.layout = FALSE, main = "Blood vessel malform.",
   xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)")

par(oldPar)


###################################################
### code chunk number 77: 07-long_prof_by_group-preparation
###################################################
ips <- list()
for (k in 1:K){
  ips[[k]] <- getProfiles(t = "month", 
     y = c("lbili", "platelet", "spiders", "jspiders", "groupHPD"), 
     id = "id", data = subset(PBC910, groupMed == k))
}
  
yvars <- c("lbili", "platelet", "jspiders")
fit.yvars <- c("lbili", "platelet", "spiders")
ylabs <- c("Log(bilirubin)", "Platelet count", 
           "Blood vessel malform. (jittered)")


###################################################
### code chunk number 78: 07-long_prof_by_group
###################################################
oldPar <- par(mfrow = c(3, 2), mar = c(5, 4, 2, 0) + 0.1, bty = "n", cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)
for (v in 1:length(yvars)){
  for (k in 1:2){
    YLIM <- range(PBC910[, yvars[v]], na.rm = TRUE)
    plotProfiles(ip = ips[[k]], 
       data = subset(PBC910, groupMed == k),
       var = yvars[v], tvar = "month", ylim = YLIM,
       gvar = "groupHPD", col = if (v <= 2) GCOL else GCOL2, bg = GCOL,
       xlab = "Time (months)",
       ylab = ifelse(k == 1, ylabs[v], ""),
       yaxt = ifelse(k == 1, "s", "n"),
       lines = (v <= 2), points = (v == 3), cex.points = 1.4,
       auto.layout = FALSE, main = "")
    lines(tpred, fit[[fit.yvars[v]]][, k], col = clCOL[k], lwd = 3)    
    if (v == 1) title(main = paste("Group", k), cex.main = 1.5)
  }  
}  
par(oldPar)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### SECTION 3.8                                                     #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 79: PED and related quantities
###################################################
print(mod$PED)


###################################################
### code chunk number 80: First ten sampled values of the model deviance - chain 1
###################################################
print(mod$Deviance1[1:10], digits = 9)


###################################################
### code chunk number 81: Running models for different numbers of mixture components (eval = FALSE)
###################################################
Devs1 <- Devs2 <- list()
PED <- data.frame()
for (K in 1:4){
  cat("Calculating K = ", K, "\n========================\n", sep="")
  if (K == 2){
    modK <- mod
  }else{    
    set.seed(20042005 + K)
    modK <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
       dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
       id = PBC910[, "id"],
       x = list(lbili    = "empty", 
               platelet = "empty", 
               spiders  = PBC910[, "month"]),
       z = list(lbili    = PBC910[, "month"], 
                platelet = PBC910[, "month"], 
                spiders  = "empty"),
       random.intercept = rep(TRUE, 3),
       prior.b = list(Kmax = K), 
       nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),
       parallel = FALSE)
  }  
  Devs1[[K]] <- modK[["Deviance1"]]
  Devs2[[K]] <- modK[["Deviance2"]]
  PED <- rbind(PED, modK[["PED"]])
  colnames(PED) <- names(modK[["PED"]])
  
  rm(list = "modK")  
}


###################################################
### code chunk number 83: Penalized expected deviances for models with different values of K
###################################################
print(PED, digits = 6)


###################################################
### code chunk number 84: Summary for the difference in deviances of models with 2 and 1 mixture components
###################################################
sD21 <- summaryDiff(c(Devs1[[2]], Devs2[[2]]), c(Devs1[[1]], Devs2[[1]]))
print(sD21, digits = 4)


###################################################
### code chunk number 85: Summary for the difference in deviances among all pairs of models
###################################################
sDiff <- data.frame()
for (K1 in 1:3){
  for (K2 in (K1 + 1):4){
    tmp1 <- summaryDiff(c(Devs1[[K2]], Devs2[[K2]]), 
                        c(Devs1[[K1]], Devs2[[K1]]))
    tmp2 <- as.data.frame(matrix(c(tmp1$Pcut, tmp1$summary), nrow = 1))
    colnames(tmp2) <- c(names(tmp1$Pcut), names(tmp1$summary))
    rownames(tmp2) <- paste(K2, "-", K1)
    sDiff <- rbind(sDiff, tmp2)
  }  
}
print(sDiff)


###################################################
### code chunk number 86: REPEAT: Print result of the previous calculation
###################################################
print(sDiff)


###################################################
### code chunk number 88: 08-cdf_deviance
###################################################
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1) + 0.1, bty = "n")
COL <- terrain_hcl(4, c = c(65, 15), l = c(45, 80), power = c(0.5, 1.5))
plot(c(14000, 14275), c(0, 1), type="n", 
     xlab="Deviance", ylab="Posterior CDF")
for (K in 1:4){
  medDEV <- median(c(Devs1[[K]], Devs2[[K]]))
  ECDF <- ecdf(c(Devs1[[K]], Devs2[[K]]))
  plot(ECDF, col = COL[K], lwd = 2, add = TRUE)
  text(medDEV + 0.5, 0.5, labels = K)
}  


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### Code from the Appendices which are available only                                     #######
####### in the extended version of the paper,                                                 #######
####### see http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf               #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### APPENDIX A                                                      #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 89: Elements of object of class GLMM_MCMClist
###################################################
names(mod)


###################################################
### code chunk number 90: Elements of object of class GLMM_MCMC
###################################################
names(mod[[1]])


###################################################
### code chunk number 91: Shift vector and scale matrix for both chains
###################################################
print(mod[[1]]$scale.b)
print(mod[[2]]$scale.b)


###################################################
### code chunk number 92: Shift vector and scale matrix for chain 1
###################################################
print(mod[[1]]$scale.b)


###################################################
### code chunk number 93: Prior for mixture related parameters for both chains
###################################################
print(mod[[1]]$prior.b)
print(mod[[2]]$prior.b)


###################################################
### code chunk number 94: Prior for fixed effects for both chains
###################################################
print(mod[[1]]$prior.alpha)
print(mod[[2]]$prior.alpha)


###################################################
### code chunk number 95: Prior for dispersion parameters for both chains
###################################################
print(mod[[1]]$prior.eps)
print(mod[[2]]$prior.eps)


###################################################
### code chunk number 96: Initial values for random effects and mixture parameters - chain 1
###################################################
print(mod[[1]]$init.b)


###################################################
### code chunk number 97: Initial values for random effects and mixture parameters
###################################################
print(mod[[2]]$init.b)


###################################################
### code chunk number 98: initial values for random effects - chain 2
###################################################
print(mod[[2]]$init.b$b[1:5,])


###################################################
### code chunk number 99: Initial values for shifted-scaled mixture means - chain 2
###################################################
print(mod[[2]]$init.b$mu)


###################################################
### code chunk number 100: Initial values for shifted-scaled mixture covariance matrices - chain 2
###################################################
print(mod[[2]]$init.b$Sigma)


###################################################
### code chunk number 101: initial values for the hyperparameter gamma.b - chain 2
###################################################
print(mod[[2]]$init.b$gammaInv)


###################################################
### code chunk number 102: initial values for fixed effects - both chains
###################################################
print(mod[[1]]$init.alpha)
print(mod[[2]]$init.alpha)


###################################################
### code chunk number 103: initial values for dispersion parameters
###################################################
print(mod[[1]]$init.eps)
print(mod[[2]]$init.eps)


###################################################
### code chunk number 104: Full specification of the prior hyperparameters by the user and restart of the MCMC simulation (eval = FALSE)
###################################################
set.seed(20072011)
modContinue <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = PBC910[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = PBC910[, "month"]),
    z = list(lbili    = PBC910[, "month"], 
             platelet = PBC910[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    scale.b = list(shift = c( 0.31516,  0.00765, 5.52621, 
                             -0.00663, -2.74954),
                   scale = c(0.8645, 0.0201, 0.3486, 0.0157, 3.2285)),                       
    prior.b = list(Kmax = 2, priormuQ = "independentC", 
                   delta = 1, xi = rep(0, 5), D = diag(rep(36, 5)), 
                   zeta = 6, gD = rep(0.2, 5), hD = rep(0.278, 5)), 
    prior.alpha = list(mean = 0, var = 10000),                       
    prior.eps = list(zeta = 2, g = 0.2, h = 2.76),
    init.b  = mod[[1]]$state.last.b,
    init2.b = mod[[2]]$state.last.b,                   
    init.alpha  = mod[[1]]$state.last.alpha,
    init2.alpha = mod[[2]]$state.last.alpha,                   
    init.eps = mod[[1]]$state.last.eps,
    init2.eps = mod[[2]]$state.last.eps,                   
    nMCMC = c(burn = 0, keep = 1000, thin = 10, info = 100),
    parallel = TRUE)


####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######
####### APPENDIX B                                                      #######
####### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #######

###################################################
### code chunk number 105: Print part of w sample
###################################################
print(mod[[1]]$w_b[1:3,])


###################################################
### code chunk number 106: Print part of shifted-scaled mu sample
###################################################
print(mod[[1]]$mu_b[1:3,])


###################################################
### code chunk number 107: Print part of scaled D sample
###################################################
print(mod[[1]]$Sigma_b[1:3,])


###################################################
### code chunk number 108: Print part of scaled Q and Li samples
###################################################
print(mod[[1]]$Q_b[1:3,])
print(mod[[1]]$Li_b[1:3,])


###################################################
### code chunk number 109: Chains with relabeled w sample
###################################################
wSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "w_b")
print(wSamp[1:3,])


###################################################
### code chunk number 110: Chains with relabeled mu sample
###################################################
muSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "mu_b")
print(muSamp[1:3,])


###################################################
### code chunk number 111: Chains with relabeled var sample
###################################################
varSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "var_b")
print(varSamp[1:3,])


###################################################
### code chunk number 112: Chains with relabeled sd sample
###################################################
sdSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "sd_b")
print(sdSamp[1:3,])


###################################################
### code chunk number 113: Chains with relabeled cor sample
###################################################
corSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "cor_b")
print(corSamp[1:3,])


###################################################
### code chunk number 114: Chains with relabeled D sample
###################################################
DSamp <- NMixChainComp(mod[[1]], relabel = TRUE, param = "Sigma_b")
print(DSamp[1:3,])


###################################################
### code chunk number 115: Print part of alpha sample
###################################################
print(mod[[1]]$alpha[1:3,])


###################################################
### code chunk number 116: Print part of sigma sample
###################################################
print(mod[[1]]$sigma_eps[1:3,])


###################################################
### code chunk number 117: Print parts of hyperparameter sample
###################################################
print(mod[[1]]$gammaInv_b[1:3,])
print(mod[[1]]$gammaInv_eps[1:3,])


###################################################
### code chunk number 118: Print part of mixture.b sample
###################################################
print(mod[[1]]$mixture_b[1:3,])


###################################################
### code chunk number 119: Print part of Deviance
###################################################
print(mod[[1]]$Deviance[1:10], digits = 6)
