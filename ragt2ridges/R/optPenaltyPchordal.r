optPenaltyPchordal <- function (Y, lambdaMin, lambdaMax, lambdaInit=(lambdaMin+lambdaMax)/2, zeros, cliques=list(), separators=list(),  target=default.target(covML(Y)), type="Alt"){ 
	##############################################################################################
	# determines the optimal value of the penalty parameter by application of the Brent 
	# algorithm to the (leave-one-out) cross-validated log-likelihood
	##############################################################################################    

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0) { stop("Matrix Y contains missings.") }
	if (as.character(class(lambdaMin)) != "numeric"){ stop("Input (lambdaMin) is of wrong class") }
	if (length(lambdaMin) != 1){ stop("lambdaMin must be a scalar") }
	if (lambdaMin <= 0){ stop("lambdaMin must be positive") }
	if (class(lambdaMax) != "numeric"){ stop("Input (lambdaMax) is of wrong class") }
	if (length(lambdaMax) != 1){ stop("lambdaMax must be a scalar") }
	if (lambdaMax <= lambdaMin){ stop("lambdaMax must be larger than lambdaMin") }
	if (as.character(class(lambdaInit)) != "numeric"){ stop("Input (lambdaInit) is of wrong class") }
	if (length(lambdaInit) != 1){ stop("lambdaInit must be a scalar") }
	if (lambdaInit <= lambdaMin){ stop("lambdaInit must be larger than lambdaMin") }
	if (lambdaMax <= lambdaInit){ stop("lambdaInit must be smaller than lambdaMax") }    
	if (!(type %in% c("Alt", "ArchI", "ArchII"))){ stop("type should be one of {'Alt', 'ArchI', 'ArchII'}") }
	if (!is.list(cliques)){ stop("Input (cliques) is of wrong class.") }    
	if (!is.list(separators)){ stop("Input (separators) is of wrong class.") }    

	# if chordal decomposition not supplied as a clique and separator list, make it so
	if (length(cliques) == 0){
		supportInfo <- support4ridgeP(nNodes=ncol(Y), zeros=zeros);
		cliques <- supportInfo$cliques; separators <- supportInfo$separators;
	}

	# determine optimal value of ridge penalty parameter
	optLambda <- optim(lambdaInit, .cvlPchordal, method="Brent", lower=lambdaMin, upper=lambdaMax, Y=Y, target=target, zeros=zeros, cliques=cliques, separators=separators, type=type)$par
	return(optLambda)
}


