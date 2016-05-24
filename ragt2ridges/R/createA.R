createA <- function(p, topology, nonzeroA=0, nCliques=1, nHubs=1, nBands=1, percZeros=0.9, stationary=TRUE){
	#####################################################################################################
	# 
	# DESCRIPTION:
	# Generates regression coefficient matrices of the VAR(1) with various type of topologies.
	#
	# ARGUMENTS:
	# -> p          : A positive \code{integer} specifying the dimension of the square matrix A.
	# -> topology   : Topology to impose on A: a \code{character} equalling either \code{"clique"},  \code{"hub"},  \code{"chain"}, or \code{"random"}.   
	# -> nonzeroA   : Numeric, value that nonzero elements of A will assume. If equal to zero, a random value from the interval [-1,1] is sampled. 
	# -> nCliques	: When \code{topology="clique"}, this positive integer specifies number of cliques.
	# -> nHubs	    : When \code{topology="hub"}, this positive integer specifies number of hubs.
	# -> nBands  	: When \code{topology="chain"}, this positive integer specifies number of bands.
	# -> percZeros	: When \code{topology="random"}, the probability with which zero elements of A are to be sampled.
	# -> stationary	: A \code{logical}: should the generated A be stationary?
	# 
	# DEPENDENCIES: 
	# require("Matrix")          # functions from package : bandSparse
	#
	#####################################################################################################

	# iput checks
	if (as.character(class(p)) != "numeric"){ stop("Input (p) is of wrong class.") }
	if (length(p) != 1){ stop("Input (p) is of wrong length.") }
	if (is.na(p)){ stop("Input (p) is not a positive integer.") }
	if (p < 0){ stop("Input (p) is not a positive integer.") }
	if (as.character(class(topology)) != "character"){ stop("Input (topology) is of wrong class.") }
	if (as.character(class(topology)) == "character"){ if (!(topology %in% c("clique", "chain", "hub", "random"))){ stop("Input (topology) ill-specified.") } }
	if (as.character(class(nonzeroA)) != "numeric"){ stop("Input (nonzeroA) is of wrong class.") }
	if (length(nonzeroA) != 1){ stop("Input (nonzeroA) is of wrong length.") }
	if (is.na(nonzeroA)){ stop("Input (nonzeroA) is not a non-negative number.") }
	if (as.character(class(nCliques)) != "numeric"){ stop("Input (nCliques) is of wrong class.") }
	if (length(nCliques) != 1){ stop("Input (nCliques) is of wrong length.") }
	if (is.na(nCliques)){ stop("Input (nCliques) is not a positive integer.") }
	if (nCliques < 0){ stop("Input (nCliques) is not a positive integer.") }
	if (nCliques > p){ stop("Input (nCliques) is not smaller than (or equal to) p.") }
	if (as.character(class(nHubs)) != "numeric"){ stop("Input (nHubs) is of wrong class.") }
	if (length(nHubs) != 1){ stop("Input (nHubs) is of wrong length.") }
	if (is.na(nHubs)){ stop("Input (nHubs) is not a positive integer.") }
	if (nHubs < 0){ stop("Input (nHubs) is not a positive integer.") }
	if (nHubs > p){ stop("Input (nHubs) is not smaller than (or equal to) p.") }
	if (as.character(class(nBands)) != "numeric"){ stop("Input (nBands) is of wrong class.") }
	if (length(nBands) != 1){ stop("Input (nBands) is of wrong length.") }
	if (is.na(nBands)){ stop("Input (nBands) is not a positive integer.") }
	if (nBands < 0){ stop("Input (nBands) is not a positive integer.") }
	if (nBands > p){ stop("Input (nBands) is not smaller than (or equal to) p.") }
	if (as.character(class(percZeros)) != "numeric"){ stop("Input (percZeros) is of wrong class.") }
	if (length(percZeros) != 1){ stop("Input (percZeros) is of wrong length.") }
	if (is.na(percZeros)){ stop("Input (percZeros) is not a non-negative number.") }
	if (percZeros <= 0){ stop("Input (percZeros) is not a positive number.") }
	if (percZeros >= 1){ stop("Input (percZeros) is not smaller than one.") }
	if (as.character(class(stationary)) != "logical"){ stop("Input (stationary) is of wrong class.") }

    # start actual function
    again <- TRUE
    while (again){

        # select a random value if the value of the nonzero elements of A has not been specified 
        if (nonzeroA==0){ nonzeroA <- runif(1, -1, 1) }

        # generate an A with a chain topology 
        if (topology == "chain") {
            diags <- list()    
            for (d in 0:nBands){
                diags[[d+1]] <- rep(nonzeroA, p-d)
            }
            A <- as.matrix(bandSparse(p, k = -c(0:nBands), diagonals=diags, symmetric=FALSE))
        }

        # generate an A with a hub topology 
        if (topology == "hub") {
            if (p %% nHubs == 0){
                hubIDs <- 1 + p / nHubs * c(0:(nHubs-1))
            } else {
                hubIDs <- 1 + floor(p / nHubs) * c(0:(nHubs-1))
            }
            A <- matrix(0, p, p)
            A[,hubIDs] <- nonzeroA
            A[upper.tri(A)] <- 0
        }

        # generate an A with a clique topology 
        if (topology == "clique") {
            if (p %% nCliques == 0){
                cliqueSizes <- rep(p / nCliques, nCliques)
            } else {
                cliqueSizes <- rep(floor(p / nCliques), nCliques)
                cliqueSizes[nCliques] <- cliqueSizes[nCliques] + p %% nCliques
            }
            A <- as.matrix(bdiag(lapply(cliqueSizes, function(x){ matrix(nonzeroA, x, x) })))
            A[upper.tri(A)] <- 0
        }

        # generate an A with a random topology 
        if (topology == "random") {
            A <- matrix(sample(0:1, p*p, replace=TRUE, prob=c(percZeros, 1-percZeros)), nrow=p, ncol=p)
            A[A!=0] <- nonzeroA
        }

        # assess stationarity   
        evs <- abs(eigen(A, only.values = TRUE)$values)
        if (max(evs) < 1 || !stationary){
            again <- FALSE
        } else {
            print("non-stationary A generated: trying again")                
        }
    } 
   
    return(A)
}


