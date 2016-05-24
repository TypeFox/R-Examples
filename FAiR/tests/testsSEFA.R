notests <- TRUE
if(notests) q(save = "no")
stopifnot(require(FAiR))

## Globals
tol <- 0.006

## Compare EFA to SEFA
manifest_Harman74 <- make_manifest(covmat = Harman74.cor)
res0_Harman74 <- make_restrictions(manifest = manifest_Harman74, factors = 3,
					model = "EFA")
EFA1_Harman74 <- Factanal(manifest = manifest_Harman74, restrictions = res0_Harman74)

x <- matrix(NA_real_, nrow(cormat(manifest_Harman74)), ncol = 3)
free <- is.na(x)
beta <- new("parameter.coef.SEFA", x = x, free = free, num_free = sum(free))

x <- diag(3)
free <- lower.tri(x)
Phi <- new("parameter.cormat", x = x, free = free, num_free = sum(free))

res1_Harman74 <- make_restrictions(manifest = manifest_Harman74, beta = beta, Phi = Phi)
SEFA1_Harman74 <- Factanal(manifest = manifest_Harman74, restrictions = res1_Harman74)

deviance(EFA1_Harman74)
deviance(SEFA1_Harman74)

RAM <- FA2RAM(SEFA1_Harman74)
library(sem)
(sem_CFA <- sem(RAM, cormat(manifest_Harman74), N = manifest_Harman74@n.obs))
deviance(SEFA1_Harman74)
df.residual(SEFA1_Harman74)
FAiR:::FAiR_stress_test(SEFA1_Harman74)

## personality example, also test DWLS
library(psych); data(bfi)
for(i in 1:ncol(bfi)) bfi[,i] <- factor(bfi[,i], ordered = TRUE)
man <- make_manifest(bfi[,1:25])
order1 <- matrix(5, 25, 5)
order2 <- matrix(25, 25, 5)

order1[1:5,1] <- order1[6:10,2] <- order1[11:15,3] <- order1[16:20,4] <-
		order1[21:25,5] <- 1

order2[1:5,1] <- order2[6:10,2] <- order2[11:15,3] <- order2[16:20,4] <-
		order2[21:25,5] <- 5

fixed1 <- matrix(0, 25, 5)
fixed1[1:5,1] <- fixed1[6:10,2] <- fixed1[11:15,3] <- fixed1[16:20,4] <- 
		fixed1[21:25,5] <- NA_real_

Phi <- diag(5)
free <- lower.tri(Phi)
Phi <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))


free <- is.na(fixed1)
beta <- new("parameter.coef", x = fixed1, free = free, num_free = sum(free))

res_bfi_CFA <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
				discrepancy = "ADF")

#set.seed(12345)
starts_CFA <- matrix(runif(1000 * 60, min = -.25, max = .25), nrow = 1000)
CFA1_bfi <- Factanal(man, res_bfi_CFA, starting.values = starts_CFA)
#library(sem)
#RAM <- FA2RAM(CFA1_bfi)
#sem_CFA <- sem(RAM, model.matrix(man, standardized = FALSE), man@n.obs)
#summary(sem_CFA)
#deviance(CFA1_bfi)
#FAiR:::FAiR_stress_test(CFA1_bfi)

#Phi <- cormat(CFA1_bfi)
#beta <- loadings(CFA1_bfi)
# start.vec <- c(Phi[lower.tri(Phi)], beta)

fixed1 <- NA * fixed1
free <- is.na(fixed1)
beta <- new("parameter.coef.SEFA", x = fixed1, free = free, num_free = sum(free))

res1_bfi <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			discrepancy = "ADF",
 			criteria = list("ranks_rows_1st", "ranks_cols_1st"),
 			methodArgs = list(row_ranks_1st = order1, col_ranks_1st = order2))
starts_SEFA <- matrix(runif(1000 * res1_bfi@nvars, min = -.25, max = .25), nrow = 1000)
starts_SEFA[1,] <- c(cormat(CFA1_bfi)[lower.tri(Phi@x)], coef(CFA1_bfi), log(CFA1_bfi@scale))
SEFA1_bfi <- Factanal(man, res1_bfi, starting.values = starts_SEFA)
summary(SEFA1_bfi)
#FAiR:::FAiR_stress_test(SEFA1_bfi)

## Compare restrictions.general with restrictions.2ndorder
fixed <- matrix(NA_real_, nrow = nrow(cormat(manifest_Harman74)), ncol = 4)
free <- is.na(fixed)
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free))

fixed <- matrix(NA_real_, nrow = 4, ncol = 1)
free <- is.na(fixed)
Delta <- new("parameter.coef", x = fixed, free = free, num_free = sum(free))

res2_Harman74 <- make_restrictions(manifest = manifest_Harman74, 
					      beta = beta, Delta = Delta)

SEFA2_Harman74 <- Factanal(manifest = manifest_Harman74, restrictions = res2_Harman74)

show(SEFA2_Harman74)
summary(SEFA2_Harman74)
RAM <- FA2RAM(SEFA2_Harman74)
sem_SEFA <- sem(RAM, model.matrix(manifest_Harman74, standardized = FALSE), 
				  manifest_Harman74@n.obs)
sem_SEFA
deviance(SEFA2_Harman74)
df.residual(SEFA2_Harman74)
FAiR:::FAiR_stress_test(SEFA2_Harman74)


fixed <- matrix(NA_real_, nrow = nrow(cormat(manifest_Harman74)), ncol = 5)
free <- is.na(fixed)
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free))


fixed <- matrix(NA_real_, nrow = 5, ncol = 2)
free <- is.na(fixed)
Delta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free))

Xi <- diag(2)
free <- lower.tri(Xi)
Xi <- new("parameter.cormat", x = Xi, free = free, num_free = sum(free))

res3_Harman74 <- make_restrictions(manifest = manifest_Harman74,
					beta = beta, Delta = Delta, Xi = Xi, 
					criteria = list("no_neg_suppressors_1st"),
					methodArgs = list(FC_threshold = -0.05))
res3_Harman74@Domains[1:11,1] <- -.8
res3_Harman74@Domains[1:11,2] <-  .8
set.seed(12345)
starts <- matrix(runif(1000 * res3_Harman74@nvars, max = 0.25), nrow = 1000)
SEFA3_Harman74 <- Factanal(manifest = manifest_Harman74, restrictions = res3_Harman74,
				starting.values = starts, boundary.enforcement = 1)
RAM <- FA2RAM(SEFA3_Harman74)
sem_SEFA <- sem(RAM, model.matrix(manifest_Harman74, standardized = FALSE), 
				  manifest_Harman74@n.obs)
#summary(sem_SEFA)
deviance(SEFA3_Harman74)
FAiR:::FAiR_stress_test(SEFA3_Harman74)
