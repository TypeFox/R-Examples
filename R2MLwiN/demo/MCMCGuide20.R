############################################################################
#     MLwiN MCMC Manual
#
# 20  Multilevel Factor Analysis Modelling . . . . . . . . . . . . . . . 303
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

# 20.1 Factor analysis modelling . . . . . . . . . . . . . . . . . . . . 303

# 20.2 MCMC algorithm . . . . . . . . . . . . . . . . . . . . . . . . . .304

# 20.3 Hungarian science exam . . . . . . . . . . . . . . . . . . . . . .304

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

# User's input if necessary

## Read hungary1 data
data(hungary1, package = "R2MLwiN")

round(colMeans(hungary1[, c("es_core", "biol_core", "phys_core")]), 4)
round(apply(hungary1[, c("es_core", "biol_core", "phys_core")], 2, sd), 4)
round(cor(hungary1[, c("es_core", "biol_core", "phys_core")]), 4)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | student), D = "Multivariate Normal", 
  data = hungary1))

covM1 <- matrix(, 6, 6)
colnames(covM1) <- rownames(covM1) <- c("cons.es_core", "cons.biol_core", "cons.biol_r3", "cons.biol_r4", "cons.phys_core", 
  "cons.phys_r2")
covM1[upper.tri(covM1, diag = TRUE)] <- mymodel@RP
# covM1[lower.tri(covM1)] <- t(covM1)[lower.tri(covM1)]
round(cov2cor(t(covM1)), 3)

# 20.4 A single factor Bayesian model . . . . . . . . . . . . . . . . . . 307

nfact <- 1
lev.fact <- 1
nfactcor <- 0
factcor <- NULL
# see FACT in the MCMC manual
loading <- matrix(c(1, 0, 0, 0, 0, 0, 1), ncol = 7, nrow = nfact, byrow = TRUE)
constr <- matrix(c(1, 0, 0, 0, 0, 0, 0), ncol = 7, nrow = nfact, byrow = TRUE)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | student), D = "Multivariate Normal", 
  estoptions = list(EstM = 1, fact = list(nfact = nfact, lev.fact = lev.fact, nfactcor = nfactcor, factcor = factcor, 
    loading = loading, constr = constr)), data = hungary1))

ranks <- rank(mymodel@fact.chains$scores)
plot(x = ranks, mymodel@fact.chains$scores, xlab = "rank", ylab = "factor scores", ylim = c(-2.35, 1.2))
abline(h = 0, lty = "dotted")

# 20.5 Adding a second factor to the model . . . . . . . . . . . . . . . .313

nfact <- 2
lev.fact <- c(1, 1)
nfactcor <- 0
factcor <- NULL
# see FACT in the MCMC manual
loading <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1), ncol = 7, nrow = nfact, byrow = TRUE)
constr <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0), ncol = 7, nrow = nfact, byrow = TRUE)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | student), D = "Multivariate Normal", 
  estoptions = list(EstM = 1, fact = list(nfact = nfact, lev.fact = lev.fact, nfactcor = nfactcor, factcor = factcor, 
    loading = loading, constr = constr)), data = hungary1))

scores <- mymodel@fact.chains$scores
plot(scores[, 2], scores[, 1], xlab = "Factor 2", ylab = "Factor 1")
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
# 20.6 Examining the chains of the loading estimates . . . . . . . . . . 317

loads <- mymodel@fact.chains$loadings

sixway(loads[, "load1_biol_core", drop = FALSE], acf.maxlag = 400, name = "load1.2")
sixway(loads[, "load2_biol_r3", drop = FALSE], acf.maxlag = 1500, name = "load2.3")

## burn-in: 5000, iterations=10,000

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | student), D = "Multivariate Normal", 
  estoptions = list(EstM = 1, fact = list(nfact = nfact, lev.fact = lev.fact, nfactcor = nfactcor, factcor = factcor, 
    loading = loading, constr = constr), mcmcMeth = list(burnin = 5000, iterations = 10000)), data = hungary1))

loads <- mymodel@fact.chains$loadings

sixway(loads[, "load2_biol_r3", drop = FALSE], acf.maxlag = 500, name = "load2.3")

# 20.7 Correlated factor . . . . . . . . . . . . . . . . . . . . . . . .319

nfact <- 2
lev.fact <- c(1, 1)
nfactcor <- 1
factcor <- c(1, 2, 0, 0)
# see FACT in the MCMC manual
loading <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1), ncol = 7, nrow = nfact, byrow = TRUE)
constr <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0), ncol = 7, nrow = nfact, byrow = TRUE)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | student), D = "Multivariate Normal", 
  estoptions = list(EstM = 1, fact = list(nfact = nfact, lev.fact = lev.fact, nfactcor = nfactcor, factcor = factcor, 
    loading = loading, constr = constr), mcmcMeth = list(burnin = 5000, iterations = 10000)), data = hungary1))

# 20.8 Multilevel factor analysis . . . . . . . . . . . . . . . . . . . 320


(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | school) + (1 | student), 
  D = "Multivariate Normal", data = hungary1))

covM1 <- matrix(, 6, 6)
colnames(covM1) <- rownames(covM1) <- c("cons.es_core", "cons.biol_core", "cons.biol_r3", "cons.biol_r4", "cons.phys_core", 
  "cons.phys_r2")
covM1[upper.tri(covM1, diag = TRUE)] <- mymodel@RP[22:42]
# covM1[lower.tri(covM1)] <- t(covM1)[lower.tri(covM1)]
round(cov2cor(t(covM1)), 3)

covM2 <- matrix(, 6, 6)
colnames(covM2) <- rownames(covM2) <- c("cons.es_core", "cons.biol_core", "cons.biol_r3", "cons.biol_r4", "cons.phys_core", 
  "cons.phys_r2")
covM2[upper.tri(covM2, diag = TRUE)] <- mymodel@RP[1:21]
# covM2[lower.tri(covM2)] <- t(covM2)[lower.tri(covM2)]
round(cov2cor(t(covM2)), 3)

# 20.9 Two level factor model . . . . . . . . . . . . . . . . . . . . . 321

nfact <- 2
lev.fact <- c(1, 2)
nfactcor <- 0
factcor <- NULL
# see FACT in the MCMC manual
loading <- matrix(c(1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1), ncol = 7, nrow = nfact, byrow = TRUE)
constr <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), ncol = 7, nrow = nfact, byrow = TRUE)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + (1 | school) + (1 | student), 
  D = "Multivariate Normal", estoptions = list(EstM = 1, fact = list(nfact = nfact, lev.fact = lev.fact, nfactcor = nfactcor, 
    factcor = factcor, loading = loading, constr = constr)), data = hungary1))

ranks <- rank(na.omit(mymodel@fact.chains$scores[, 2]))
plot(x = ranks, y = na.omit(mymodel@fact.chains$scores[, 2]), xlab = "rank", ylab = "factor scores", ylim = c(-1, 
  1))
abline(h = 0, lty = "dotted")

# 20.10 Extensions and some warnings . . . . . . . . . . . . . . . . . . 324

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .325





############################################################################
