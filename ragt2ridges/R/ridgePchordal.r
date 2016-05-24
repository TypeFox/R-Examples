ridgePchordal <- function(S, lambda, zeros, cliques=list(), separators=list(), target=default.target(S), type="Alt", optimizer="nlm", grad=FALSE, verbose=TRUE, ...){
    #####################################################################################################
	# 
	# DESCRIPTION:
	# Ridge estimation of the precision matrix with known zeros.
	# Nonzeros should form a chordal graph. If not chordal, the function makes it so.
	#
	# ARGUMENTS:
	# -> S          	: Sample covariance matrix.
	# -> lambda     	: A numeric representing the value of the penalty parameter.
	# -> target     	: A target matrix (in precision terms) for Type I ridge estimators.
	# -> zeros	        : A matrix with indices of entries of the precision matrix 
	#                     that are constrained to zero. The matrix comprises two columns, each 
	#                     row corresponding to an entry of the adjacency matrix. The first 
	#                     column contains the row indices and the second the column indices. The 
	#                     specified graph should be undirected and decomposable. If not, it is 
	#                     symmetrized and triangulated. Hence, it may differ from the input 'zeros'.
	# -> cliques        : A 'list'-object containing the node indices per clique as 
	#                     obtained from the 'support4ridgeP'-function.
	# -> separators     : A 'list'-object containing the node indices per separator as 
	# 				      obtained from the 'support4ridgeP'-function.
	# -> type       	: A character indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
	# -> optimizer     	: Which optimization function should be used: "optim" or "nlm"?
	# -> grad       	: Logical indicator: should, next to the precision matrix estimate, also the gradient be returned? 
	# -> verbose    	: Logical indicator: should intermediate output be printed on the screen? 
	# -> ...        	: Additional arguments passed on to either "optim" or "nlm".
	#
	# DEPENDENCIES: 
	# require("igraph")          # functions from package : graph.adjancency, igraph.to.graphNEL
	# require("gRbase")          # functions from package : triangulate
	# require("RGBL")            # functions from package : is.triangulated 
	# require("rags2ridges")     # functions from package : adjacentMat, default.target
	# require("stats")           # functions from package : nlm, optim
	# require("utils")           # functions from package : txtProgressBar, setTxtProgressBar, Sys.sleep
	#    
	# NOTES:
	# 1) Currently, uses the full inverse instead of the partial inverse. To be fixed in the far future.
	#
	#####################################################################################################

   	# input checks
	if (!is.matrix(S)){ stop("S should be a matrix") }
	if (!isSymmetric(S)){ stop("S should be a symmetric matrix") }
	if (length(lambda) != 1){ stop("Input (lambda) is of wrong length.") }
	if (is.na(lambda)){ stop("Input (lambda) is not a positive number.") }
	if (lambda <= 0){ stop("Input (lambda) is not a positive number.") }
	if (!is.null(zeros) & as.character(class(zeros)) != "matrix"){ stop("Input (zeros) is of wrong class.") }    
	if (!is.null(zeros)){ if(ncol(zeros) != 2){ stop("Wrong dimensions of the (zeros).") } } 
	if (!is.null(zeros)){ zeros <- zeros[order(zeros[,2], zeros[,1]),] }
	if (!(type %in% c("Alt", "ArchI", "ArchII"))) {  stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")   }
	if (!(optimizer %in% c("nlm", "optim"))) {  stop("type should be one of {'nlm', 'optim'}")   }
	if (as.character(class(verbose)) != "logical"){ stop("Input (verbose) is of wrong class.") }
	if (!is.list(cliques)){ stop("Input (cliques) is of wrong class.") }    
	if (!is.list(separators)){ stop("Input (separators) is of wrong class.") }    

	# intermediate output
	if(verbose){ 
		cat(" ", "\n") 
		cat("Progress report ....", "\n")
		cat(paste("-> ----------------------------------------------------------------------", sep=""), "\n") 
	}

	# if chordal decomposition not supplied as a clique and separator list, make it so
	if (length(cliques) == 0){
		supportInfo <- support4ridgeP(nNodes=nrow(S), zeros=zeros);
		cliques <- supportInfo$cliques; separators <- supportInfo$separators;
	}

	# make adjacency matrix of support
	Pinit <- matrix(1, ncol=ncol(S), nrow=nrow(S));
	Pinit[zeros] <- 0; 
	diag(Pinit) <- 1;

	# ensure the target has same support as chordal graph
	target[zeros] <- 0;
	
	# obtain the graph components
	diag(Pinit) <- 0 
	components <- clusters(graph.adjacency(Pinit, mode="undirected"))$membership
    
	# construct init estimate
	Pinit <- .armaRidgePchordalInit(S=S, lambda=lambda, target=target, type=type, Cliques=cliques, Separators=separators)
	evs <- eigen(Pinit, only.values=TRUE)$values

	if (type=="Alt"){ 
		if (any(eigen(Pinit, only.values=TRUE)$values < 0)){ Pinit <- (Pinit + diag(diag(Pinit)))/2 }	

		# number of non-converged components
		nonconvergedComp <- 0

		# report number of to be estimated parameters
		diag(Pinit) <- diag(Pinit) / 2
		Pinit[upper.tri(Pinit)] <- 0
		nonzeros <- which(Pinit != 0, arr.ind=TRUE)
		Xinit <- Pinit[nonzeros]
		if (verbose){ 
			cat(paste("-> optimization over                   : ", nrow(nonzeros), " out of ", ncol(S) * (ncol(S) + 1) / 2, " unique", sep=""), "\n") 
			cat(paste("->                                       ", "parameters (", round(200 * nrow(nonzeros) /  (ncol(S) * (ncol(S) + 1)), digits=2), "%)", sep=""), "\n") 
			cat(paste("-> cond. number of initial estimate    : ", round(max(evs) / max(min(evs), 0), digits=2), " (if >> 100, consider", sep=""), "\n")
			cat(paste("->                                       larger values of lambda)"), "\n")
			cat(paste("-> # graph components                  : ", length(unique(components)), " (optimization per component)", sep=""), "\n")
			cat(paste("-> optimization per component          : ", sep=""), "\n")  
		}	

		if (verbose){ pBar <- txtProgressBar(min=0, max=abs(length(unique(components))), style=3, char=".") }

		# estimate per graph component
		for (subG in unique(components)){
			if (verbose){ setTxtProgressBar(pBar, subG); Sys.sleep(10^(-10)) }
	
			# construct component data    
			ids <- which(components == subG)   
			Psub <- Pinit[ids, ids, drop=FALSE] 
			nonzeros <- which(Psub != 0, arr.ind=TRUE)
			x0 <- Psub[nonzeros]
			rm(Psub)

			# matrices for alternative parametrization of the sparse precision matrix
			E1 <- matrix(0, ncol=nrow(nonzeros), nrow=length(ids))
			E1[cbind(nonzeros[,1], c(1:nrow(nonzeros)))] <- 1
			E2 <- matrix(0, ncol=nrow(nonzeros), nrow=length(ids))
			E2[cbind(nonzeros[,2], c(1:nrow(nonzeros)))] <- 1

			if (nrow(nonzeros) != length(ids) * (length(ids) + 1) / 2){	        
				# minimize minus penalized log-likelihood
                if (optimizer=="optim"){
                    xhat <- optim(x0, .armaPenLLreparP, gr=.armaPenLLreparPgrad, method="BFGS", E1=E1, E2=E2, S=S[ids,ids,drop=FALSE], 
                                    lambda=lambda, target=target[ids,ids,drop=FALSE], nonzerosR=nonzeros[,1], nonzerosC=nonzeros[,2], ...) 
                    if(xhat$convergence != 0){ nonconvergedComp <- nonconvergedComp + 1 }
                    Pinit[ids, ids] <- E1 %*% diag(xhat$par, ncol=length(x0)) %*% t(E2) + E2 %*% diag(xhat$par, ncol=length(x0)) %*% t(E1)                  	    	        
                } 
	        
                if (optimizer=="nlm"){
                    xhat <- nlm(.armaPenLLreparPforNLM, p=x0, E1=E1, E2=E2, S=S[ids,ids,drop=FALSE], lambda=lambda, 
                                    target=target[ids,ids,drop=FALSE], nonzerosR=nonzeros[,1], nonzerosC=nonzeros[,2], check.analyticals=FALSE, ...)
                    if(xhat$code != 0 & xhat$code != 1){ nonconvergedComp <- nonconvergedComp + 1 }
                    Pinit[ids, ids] <- E1 %*% diag(xhat$estimate, ncol=length(x0)) %*% t(E2) + E2 %*% diag(xhat$estimate, ncol=length(x0)) %*% t(E1)
				}
            } else {           		
                Pinit[ids, ids] <- E1 %*% diag(x0, ncol=length(x0)) %*% t(E2) + E2 %*% diag(x0, ncol=length(x0)) %*% t(E1)	    			        	    
	    	}
		}
	}

	if (verbose){ 
		if (type=="Alt"){ cat("\n"); cat("\n") } 
		cat(paste("-> estimation done ...", sep=""), "\n")
		cat(paste("-> formatting output ...", sep=""), "\n")                              
	}

	# reformatting for reporting the gradient
	diag(Pinit) <- diag(Pinit) / 2
	Pinit[upper.tri(Pinit)] <- 0
	nonzeros <- which(Pinit != 0, arr.ind=TRUE)
	x0 <- Pinit[nonzeros]
	E1 <- matrix(0, ncol=nrow(nonzeros), nrow=nrow(Pinit))
	E1[cbind(nonzeros[,1], c(1:nrow(nonzeros)))] <- 1
	E2 <- matrix(0, ncol=nrow(nonzeros), nrow=nrow(Pinit))
	E2[cbind(nonzeros[,2], c(1:nrow(nonzeros)))] <- 1
	rm(Pinit)

	if (verbose & type=="Alt"){ 
		cat(paste("-> overall summary ....", sep=""), "\n") 
		cat(paste("-> initial pen. log-likelihood         : ", round(- .armaPenLLreparP(Xinit, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]), 8), sep=""), "\n") 
		cat(paste("-> optimized pen. log-likelihood       : ", round(- .armaPenLLreparP(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]), 8), sep=""), "\n")
		if (nonconvergedComp == 0){ 
			cat("-> optimization                        : converged (most likely)", "\n") 
			cat("->                                       for all components", "\n")
		} 
		if (nonconvergedComp > 0){ 
			cat("-> optimization                        : for ", nonconvergedComp, " components", "\n")
			cat("->                                       max. no. iterations reached (or", "\n")
			cat("->                                       other indication of possible", "\n")
			cat("->                                       convergence failure)", "\n")    
		}
	}
	if (verbose){ cat(paste("-> ----------------------------------------------------------------------", sep=""), "\n")  }

	# return precision matrix from alternative parametrization
	if (!grad){
        return(E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1))
	} else {		        
		if (type=="Alt"){ return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1), grad=.armaPenLLreparPgrad(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]))) }
		if (type=="ArchI"){ return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1), grad=.armaPenLLreparGradArchI(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]))) }
		if (type=="ArchII"){ return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1), grad=.armaPenLLreparGradArchII(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]))) }
	}
}

