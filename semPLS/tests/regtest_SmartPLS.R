library(semPLS)


### Data: mobi; Model: ECSImobi
data(ECSImobi)


# scheme: centroid
# Results match with those from SmartPLS (Version: 2.0.M3)
ecsi <- sempls(ECSImobi, mobi, wscheme="centroid", sum1=TRUE)

# Import SmartPLS results
ecsiSmart <- list()

ecsiSmart$outerW <-
read.table("./SmartPLS/results/ecsi/centroid/outerWeights.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$pC <-
read.table("./SmartPLS/results/ecsi/centroid/pathCoefficients.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$tpC <-
read.table("./SmartPLS/results/ecsi/centroid/totalEffects.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$crossL <-
read.table("./SmartPLS/results/ecsi/centroid/crossLoadings.txt", header=TRUE, row.names=1, dec=",")

ecsiSmart <- lapply(ecsiSmart, as.matrix)

ecsiSmart$outerW <- ecsiSmart$outerW[rownames(ecsi$outer_weights), colnames(ecsi$outer_weights)]
ecsiSmart$pC <- ecsiSmart$pC[rownames(ecsi$path_coefficients), colnames(ecsi$path_coefficients)]
ecsiSmart$total_effects <- ecsiSmart$tpC[rownames(ecsi$total_effects), colnames(ecsi$total_effects)]
ecsiSmart$crossL <- ecsiSmart$crossL[rownames(ecsi$cross_loadings), colnames(ecsi$cross_loadings)]

# compare the results

# outer weights
stopifnot(all.equal(ecsi$outer_weights, apply(ecsiSmart$outerW, 2, semPLS:::sum1)))
# path coefficients
stopifnot(all.equal(ecsi$path_coefficients, ecsiSmart$pC))
# total effects
stopifnot(all.equal(ecsi$total_effects, ecsiSmart$total_effects))
# cross loadings
stopifnot(all.equal(ecsi$cross_loadings, ecsiSmart$crossL))



# scheme: factorial
# Results match with those from SmartPLS (Version: 2.0.M3)
ecsi <- sempls(ECSImobi, mobi, wscheme="factorial", sum1=TRUE)

# Import SmartPLS results
ecsiSmart <- list()

ecsiSmart$outerW <-
read.table("./SmartPLS/results/ecsi/factorial/outerWeights.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$pC <-
read.table("./SmartPLS/results/ecsi/factorial/pathCoefficients.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$tpC <-
read.table("./SmartPLS/results/ecsi/factorial/totalEffects.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$crossL <-
read.table("./SmartPLS/results/ecsi/factorial/crossLoadings.txt", header=TRUE, row.names=1, dec=",")

ecsiSmart <- lapply(ecsiSmart, as.matrix)

ecsiSmart$outerW <- ecsiSmart$outerW[rownames(ecsi$outer_weights), colnames(ecsi$outer_weights)]
ecsiSmart$pC <- ecsiSmart$pC[rownames(ecsi$path_coefficients), colnames(ecsi$path_coefficients)]
ecsiSmart$total_effects <- ecsiSmart$tpC[rownames(ecsi$total_effects), colnames(ecsi$total_effects)]
ecsiSmart$crossL <- ecsiSmart$crossL[rownames(ecsi$cross_loadings), colnames(ecsi$cross_loadings)]

# compare the results
# outer weights
stopifnot(all.equal(ecsi$outer_weights, apply(ecsiSmart$outerW, 2, semPLS:::sum1)))
# path coefficients
stopifnot(all.equal(ecsi$path_coefficients, ecsiSmart$pC))
# total effects
stopifnot(all.equal(ecsi$total_effects, ecsiSmart$total_effects))
# cross loadings
stopifnot(all.equal(ecsi$cross_loadings, ecsiSmart$crossL))



# scheme: path weighting
# Results match with those from SmartPLS (Version: 2.0.M3)
ecsi <- sempls(ECSImobi, mobi, wscheme="pathWeighting", sum1=TRUE)

# Import SmartPLS results
ecsiSmart <- list()

ecsiSmart$outerW <-
read.table("./SmartPLS/results/ecsi/pathWeighting/outerWeights.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$pC <-
read.table("./SmartPLS/results/ecsi/pathWeighting/pathCoefficients.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$tpC <-
read.table("./SmartPLS/results/ecsi/pathWeighting/totalEffects.txt", header=TRUE, row.names=1, dec=",")
ecsiSmart$crossL <-
read.table("./SmartPLS/results/ecsi/pathWeighting/crossLoadings.txt", header=TRUE, row.names=1, dec=",")

ecsiSmart <- lapply(ecsiSmart, as.matrix)

ecsiSmart$outerW <- ecsiSmart$outerW[rownames(ecsi$outer_weights), colnames(ecsi$outer_weights)]
ecsiSmart$pC <- ecsiSmart$pC[rownames(ecsi$path_coefficients), colnames(ecsi$path_coefficients)]
ecsiSmart$total_effects <- ecsiSmart$tpC[rownames(ecsi$total_effects), colnames(ecsi$total_effects)]
ecsiSmart$crossL <- ecsiSmart$crossL[rownames(ecsi$cross_loadings), colnames(ecsi$cross_loadings)]

# compare the results
# outer weights
stopifnot(all.equal(ecsi$outer_weights, apply(ecsiSmart$outerW, 2, semPLS:::sum1)))
# path coefficients
stopifnot(all.equal(ecsi$path_coefficients, ecsiSmart$pC))
# total effects
stopifnot(all.equal(ecsi$total_effects, ecsiSmart$total_effects))
# cross loadings
stopifnot(all.equal(ecsi$cross_loadings, ecsiSmart$crossL))
