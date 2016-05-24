CIGofVAR1 <- function(sparseA, sparseP, type="global"){
    
	############################################################################################################
	#
	# DESCRIPTION:
	# Constructs the global or contemporaneous conditional independence graph (CIG) of the VAR(1) model.
	#
	# ARGUMENTS:
	# -> sparseA                 : Matrix A of regression parameters, which is assumed to be sparse.
	# -> sparseP                 : Matrix P of precision of the error, which is assumed to be sparse.
	# -> type                    : A 'character' indicating whether the 'global' or 'contemp' (contemporaneous)
	#                              adjanceny matrix of the conditional independence graph is requested.
	# 
	# DEPENDENCIES:
	# ....
	#
	# NOTES:
	# ....
	#
	# REFERENCES: 
	# -> Dahlhaus (2000), "Graphical interaction models for multivariate time series", Metrika, 51, 157-172.
	# -> Dahlhaus, Eichler (2003), "Causality and graphical models in time series analysis", Oxford Statistical Science Series, 115-137. 
	#     	
	############################################################################################################

	# input checks    
  	if (as.character(class(sparseA)) != "matrix"){ stop("Input (sparseA) is of wrong class.") }
	if (nrow(sparseA) != ncol(sparseA)){ stop("Matrix sparseA is not square.") }
	if (as.character(class(sparseP)) != "matrix"){ stop("Input (sparseP) is of wrong class.") }
	if (nrow(sparseP) != ncol(sparseP)){ stop("Matrix sparseP is not square.") }
	if (nrow(sparseA) != ncol(sparseP)){ stop("Matrix sparseA and sparseP are of equal dimensions.") }
  	if (as.character(class(type)) != "character"){ stop("Input (type) is of wrong class.") }
  	if (!(type %in% c("global", "contemp"))){ stop("Input (type) ill-specified.") }

	# construct adjacency matrix of global markov (in)dependencies
	if (type=="global"){
		CIG <- abs(t(sparseA) %*% sparseP) + abs(sparseP %*% sparseA) + abs(t(sparseA) %*% sparseP %*% sparseA)
		diag(sparseP) <- 0
		CIG <- CIG + abs(sparseP) 
		CIG[CIG != 0] <- 1
	}
    
	# construct adjacency matrix of global markov (in)dependencies
	if (type=="contemp"){
		diag(sparseP) <- 0
		CIG <- abs(sparseP) 
		CIG[CIG != 0] <- 1
	}
	return(CIG)
}

