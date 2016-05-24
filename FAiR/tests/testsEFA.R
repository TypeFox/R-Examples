notests <- TRUE
if(notests) q(save = "no")
stopifnot(require(FAiR))

## Globals
tol <- 0.006

## Compare factanal() with Factanal()
efa_Harman23 <- factanal(covmat = Harman23.cor, factors = 2, rotation = "none")

manifest_Harman23 <- make_manifest(covmat = Harman23.cor)
res_Harman23 <- make_restrictions(manifest = manifest_Harman23, factors = 2, model="EFA")

EFA_Harman23 <- Factanal(manifest = manifest_Harman23, restrictions = res_Harman23)

show(EFA_Harman23)
summary(EFA_Harman23)

stopifnot(all.equal(efa_Harman23$uniquenesses, uniquenesses(EFA_Harman23), tol = tol))
RAM <- FA2RAM(EFA_Harman23)
library(sem)
RAM
sem_Harman23 <- sem(RAM, cormat(manifest_Harman23), manifest_Harman23@n.obs)
sem_Harman23
deviance(EFA_Harman23)
df.residual(EFA_Harman23)
FAiR:::FAiR_stress_test(EFA_Harman23)

EFA_Harman23 <- Factanal(manifest = manifest_Harman23, restrictions = res_Harman23,
			impatient = TRUE)
stopifnot(all.equal(efa_Harman23$uniquenesses, uniquenesses(EFA_Harman23), tol = tol))
FAiR:::FAiR_stress_test(EFA_Harman23)

## Check Rotation() criteria
EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("geomin"),
			methodArgs = list(nfc_threshold = 0.25, delta = .01, 
					matrix = "PP"),	max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("quartimin"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "FC"), 
				max.generations = 100, solution.tolerance = 1e-4)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("target"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "PP",
						Target = matrix(c(      1, 0,
									1, 0,
									1, 0,
									1, 0,
									0, 1,
									0, 1,
									0, 1,
									1/4, 3/4), 
								ncol = 2, byrow = TRUE)), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("pst"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "PP",
						Target = matrix(c(      NA, 0,
									NA, 0,
									NA, 0,
									NA, 0,
									0, NA,
									0, NA,
									0, NA,
									NA, NA), 
								ncol = 2, byrow = TRUE)), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("oblimax"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "PP"), 
				max.generations = 100, solution.tolerance = 8e-2)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("simplimax"),
				methodArgs = list(nfc_threshold = 0.25, k = 8,
						matrix = "PP"), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("bentler"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "RS"), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("cf"),
				methodArgs = list(nfc_threshold = 0.25, kappa = 0,
						matrix = "PP"), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("infomax"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "PP"), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("mccammon"),
				methodArgs = list(nfc_threshold = 0.25, matrix = "PP"), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("oblimin"),
				methodArgs = list(nfc_threshold = 0.25, gam = 0, 
						matrix = "PP"), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("minimaximin"),
				methodArgs = list(nfc_threshold = 0.25),
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("LS"),
				methodArgs = list(nfc_threshold = 0.25),
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

EFA_Harman23_rotated <- Rotate(EFA_Harman23, criteria = list("phi"),
				methodArgs = list(nfc_threshold = 0.25, c = 1),
				max.generations = 100)
show(EFA_Harman23_rotated)
summary(EFA_Harman23_rotated)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

## Check that constraints work
# limit primary factor correlations
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"limit_correlations", "phi"),
				methodArgs = list(nfc_threshold = 0.25, lower = -.9,
						upper = .9, c = 1), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# positive manifold (without restricting factor correlations)
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"positive_manifold", "phi"),
				methodArgs = list(nfc_threshold = 0.25, 
						pm_threshold = -.1, c = 1), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# ranks_rows_1st
ranks <- matrix(2, nrow = 8, ncol = 2)
ranks[1,1] <- 1
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"ranks_rows_1st", "phi"),
				methodArgs = list(nfc_threshold = 0.25, row_ranks = ranks,
						c = 1), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# ranks_cols_1st
ranks <- matrix(8, nrow = 8, ncol = 2)
ranks[1,1] <- 7
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"ranks_cols_1st", "phi"),
				methodArgs = list(nfc_threshold = 0.25, col_ranks = ranks,
						c = 1), max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# indicators_1st
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"indicators_1st", "phi"),
				methodArgs = list(nfc_threshold = 0.25, c = 1,
						indicators = c(2,5)), 
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# evRF_1st
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"evRF_1st", "phi"),
				methodArgs = list(nfc_threshold = 0.25, c = 1),
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

# evPF_1st
EFA_Harman23_rotated <- Rotate(EFA_Harman23, 
				criteria = list("no_factor_collapse",
						"evPF_1st", "phi"),
				methodArgs = list(nfc_threshold = 0.25, c = 1),
				max.generations = 100)
FAiR:::FAiR_stress_test(EFA_Harman23_rotated)

## Compare restrictions.factanal() with restrictions.orthonormal()
efa_ability.cov <- factanal(covmat = ability.cov, factors = 2, rotation = "none")
manifest_ability.cov <- make_manifest(covmat = ability.cov)
x <- matrix(NA_real_, 6, 2)
rownames(x) <- rownames(cormat(manifest_ability.cov))
x[1,2] <- 0
free <- is.na(x)
beta <- new("parameter.coef", x = x, free = free, num_free = sum(free))
res_ability.cov <- make_restrictions(manifest = manifest_ability.cov, 
					beta = beta, discrepancy = "MLE")

EFA_ability.cov <- Factanal(manifest = manifest_ability.cov,
			restrictions = res_ability.cov)
show(EFA_ability.cov)
summary(EFA_ability.cov)
stopifnot(all.equal(efa_ability.cov$uniquenesses, 
                    EFA_ability.cov@uniquenesses, tol = tol))
FAiR:::FAiR_stress_test(EFA_ability.cov)
