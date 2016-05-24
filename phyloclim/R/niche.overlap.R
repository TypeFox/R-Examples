niche.overlap <- function(x){
	
	# CASE 1: x is a pno matrix
	# -------------------------
	if ( is.data.frame(x) ){
		x <- x[, -1]
		nspec <- ncol(x)
		DI <- matrix(nrow = nspec, ncol = nspec)
		rownames(DI) <- colnames(DI) <- names(x)
		for (i in 1:(nspec - 1)){
			for (j in (i + 1):nspec){
				dhi <- di.pno(x = x[, i], y = x[, j])
				DI[i, j] <- dhi["D"]
				DI[j, i] <- dhi["I"]
			}
		}
	}

	## CASE 2: is vector of filenames
	## ------------------------------
	if ( class(x) == "character" ){
		nspec <- length(x)
		DI <- matrix(nrow = nspec, ncol = nspec)
		rownames(DI) <- colnames(DI) <- x
		for (i in 1:(nspec - 1)){
			X <- read.asciigrid(x[i])
			for (j in (i + 1):nspec){
				Y <- read.asciigrid(x[j])
				dhi <- di.enm(x = X, y = Y)
				DI[i, j] <- dhi[1]
				DI[j, i] <- dhi[2]
			}
		}
	}
	
	## CASE 3: x is a list of 'SpatialGrid' objects
	## --------------------------------------------
	if (inherits(x[[1]], "SpatialGrid")){
		nspec <- length(x)
		DI <- matrix(nrow = nspec, ncol = nspec)
		rownames(DI) <- colnames(DI) <- names(x)
    system.time(
		for (i in 1:(nspec - 1)){
			for (j in (i + 1):nspec){
				dhi <- di.enm(x = x[[i]], y = x[[j]])
				DI[i, j] <- dhi[1]
				DI[j, i] <- dhi[2]
			}
		}
		) # system.time
	}
	class(DI) <- "niolap"
	DI
}