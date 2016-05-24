## Define class unions for optional slots, e.g. for definition
##  of slots which will be computed on demand, like the
##  mahalanobis/robust distances
## 
## "Uvector", "Umatrix", "Ulist", "Ufunction", "Utable", "UCovControl"
## are already defined in rrcov-AllClasses.R
setClassUnion("Unumeric", c("numeric", "NULL"))
setClassUnion("Ulogical", c("logical", "NULL"))


###################### FA ####################################
setClass("Fa", representation(call="language",
			      converged="Ulogical",
			      loadings="matrix",
			      communality="Uvector",
			      uniquenesses="vector",
				  cor="Ulogical",
				  covariance="matrix",
			      correlation="matrix",
				  usedMatrix="matrix",
				  reducedCorrelation="Umatrix",
			      criteria="Unumeric",
			      factors="numeric",
			      dof="Unumeric",
			      method="character",
			      scores="Umatrix",
			      scoresMethod="character",
			      scoringCoef="Umatrix",
			      meanF="Uvector",
			      corF="Umatrix",
			      STATISTIC="Unumeric",
			      PVAL="Unumeric",
			      n.obs="numeric",
			      center="Uvector",
			      eigenvalues="vector",
			      cov.control="UCovControl",
                              "VIRTUAL")) # VIRTUAL class

setClass("SummaryFa", representation(faobj = "Fa",
                                      importance  ="matrix"))

setClass("FaClassic", contains="Fa")

setClass("FaRobust", representation("VIRTUAL"),
                    contains="Fa") # VIRTUAL class

setClass("FaCov", representation(),
                    contains="FaRobust")