################################
##	Code for examples
require(ClustMMDD)

##	Main package examples

data(genotype1)
?genotype1

data(genotype2)
?genotype2

data(genotype2_ExploredModels)
head(genotype2_ExploredModels)
?genotype2_ExploredModels

#Calibration of the penalty function
outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5, header = TRUE)
cte1 = outDimJump[[1]][1]
outSlection = model.selection.R(genotype2_ExploredModels, cte = cte1, header = TRUE)
outSlection

##	== method for modelKS class
showClass("modelKS")
slotNames("modelKS")
data(exModelKS)
exModelKS
exModelKS == exModelKS

##	Backward-explorer
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
head(genotype2)
########
# The following command create a file "genotype2_ExploredModels.txt" that containts the explored models.
#output = backward.explorer(genotype2, Kmax = 10, ploidy = 2, Kmin = 1, Criterion = "CteDim")
data(genotype2_ExploredModels)
head(genotype2_ExploredModels)

##	cutEachCol
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[,-11], ploidy = 2) #ploidy = 2 for these genotype data
head(genotype2)

##	dataR2C
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[, -11], 2) #ploidy = 2 for these genotype data
head(genotype2)
genotype3 = dataR2C(genotype2, ploidy = 2)
head (genotype3$data)
str(genotype3)

##	dimJump
data(genotype2_ExploredModels)
outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5, header = TRUE)
outDimJump[[1]]

## em.cluster.R
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
head(genotype2)
#See the EM options
EmOptions() # Options can be set by \code{\link{setEmOptions()}}
par5 = em.cluster.R (genotype2, K = 5, S = c(rep(TRUE, 8), rep(FALSE, 2)), ploidy = 2)
slotNames(par5)
head(par5["membershipProba"])
par5["mixingProportions"]
par5

##	EmOptions
EmOptions()
setEmOptions(list(epsi = 1e-6))
EmOptions()
setEmOptions() # To set default values
EmOptions()

##	exModelKS
data(exModelKS)
slotNames("modelKS")
head(exModelKS["membershipProba"])
exModelKS["mixingProportions"]
exModelKS

##	genotype1
data(genotype1)
head(genotype1)

##	genotype2
data(genotype2)
head(genotype2)
data(genotype1)
genotype3 = cutEachCol(genotype1[,-11], ploidy = 2)
head(genotype3)

##	genotype2_ExploredModels
data(genotype2_ExploredModels)
head(genotype2_ExploredModels)
plot(genotype2_ExploredModels[, c("dim", "logLik")], col = "blue", xlab = "Dimension", ylab = "Log-likelihood")
# Data-driven calibration of the penalty
dimJump.R(genotype2_ExploredModels, h = 5, N=1000, header=T)

##	is.element
data(exModelKS)
is.element(c(exModelKS), c(exModelKS))
is.element(c(exModelKS, 1, c(1:5)), c(exModelKS))
is.element(c(exModelKS), c(exModelKS, 1, list(1:5, 0)))
is.element(exModelKS, c(exModelKS, 1, list(1:5, 0))) # Error

##	is.modelKS
data(exModelKS)
is.modelKS(exModelKS)
is.modelKS(1:7)

##	isInFile.R
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
head(genotype2)
	
S = c(rep(TRUE, 8), rep(FALSE, 2))
outPut = selectK.R(genotype2, S, Kmax = 6, ploidy = 2, Kmin=1)
isInFile.R(K = 5, S, "genotype2_ExploredModels.txt", header = TRUE)
isInFile.R(K = 5, rep(TRUE, 10), "genotype2_ExploredModels.txt", header = TRUE)

##	model.selection.R
data(genotype2_ExploredModels)
outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5, header = TRUE)
cte1 = outDimJump[[1]][1]
outSlection = model.selection.R(genotype2_ExploredModels, cte = cte1, header = TRUE)
outSlection

##	modelKS
data(exModelKS)
showClass("modelKS")
slotNames("modelKS")
exModelKS
exModelKS["K"]
exModelKS["S"]

##	selectK.R
data(genotype1)
head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
head(genotype2)
S = c(rep(TRUE, 8), rep(FALSE, 2))
outPut = selectK.R(genotype2, S, Kmax = 6, ploidy = 2, Kmin=1)
outPut[["BIC"]]

##	setOptions
EmOptions()
setEmOptions(list(epsi = 1e-6))
EmOptions()
setEmOptions() # To set default values
EmOptions()

##	simulData
data(exModelKS)
exModelKS
exData = simulData(exModelKS, 1000, 2)
str(exData)
head(exData$data)
head(exData$class)


## []
data(exModelKS)
slotNames(exModelKS)
exModelKS["K"]
exModelKS["S"]
exModelKS["logLik"]

