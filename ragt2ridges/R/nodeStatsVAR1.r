nodeStatsVAR1 <- function(sparseA, sparseP, as.table = FALSE){
	#####################################################################################################
	#
	# DESCRIPTION:
	# -> Function that calculates various network statistics from a sparse VAR(1) model.
	#
	# ARGUMENTS:
	# -> sparseA  : sparse regression coefficient matrix
	# -> sparseP  : sparse precision/partial correlation matrix
	# -> as.table : logical indicating if output should be returned as table; default = FALSE.
	# 
	# NOTES (network statistics produced):
	# -> degreeAin         	 : number of (temporal) edges pointing to each node ('in'-degree).
	# -> degreeAout        	 : number of (temporal) edges leaving each node ('out'-degree).
	# -> nNegAin           	 : number of negative (temporal) edges pointing to each node.
	# -> nPosAin           	 : number of positive (temporal) edges pointing to each node ('in'-degree)
	# -> nNegAout          	 : number of negative (temporal) edges leaving each node ('out'-degree)
	# -> nPosAout          	 : number of positive (temporal) edges leaving each node ('out'-degree)
	# -> degreePe          	 : number of contemporaneous edges of each node (as implied by the error precision matrix)
	# -> betweennessPe     	 : vector representing the contemporaneous betweenness centrality for each node.
	# -> closenessPe       	 : vector representing the contemporaneous closeness centrality for each node.
	# -> eigenCentralityPe 	 : vector representing the contemporaneous eigen centrality for each node.
	# -> nNegPe            	 : vector representing the number of negative contemporaneous edges for each node.
	# -> nPosPe            	 : vector representing the number of positive contemporaneous edges for each node.
	# -> variancePe        	 : vector representing the error variance of each node.
	# -> partialVarPe      	 : vector representing the partial error variance of each node.
	# -> varianceY		 : vector representing the variance of each node.
	# -> degreePy		 : number of edges of each node in the global Markov graph.
	# -> betweennessPy	 : vector representing the betweenness centrality for each node in the global Markov graph.
	# -> closenessPy	 : vector representing the closeness centrality for each node in the global Markov graph.
	# -> eigenCentralityPy	 : vector representing the eigen centrality for each node in the global Markov graph.
	# -> mutualInfo_Tplus1	 : vector with for each node its mutual information with all other nodes at the next (t+1) time point.
	# -> mutualInfo_Tplus2	 : vector with for each node its mutual information with all other nodes at the (t+2)-th time point.
	# -> itemResponse_Tplus1 : vector with for each node its  mean absolute impulse response on all other nodes at the next (t+1) time point.
	# -> itemResponse_Tplus2 : vector with for each node its  mean absolute impulse response on all other nodes at the (t+2)-th time point.
	# - Future versions of this function may include additional statistics
	# 
	# DEPENDENCIES:
	# require("igraph")          # functions from package : graph.adjacency, degree, closeness, betweenness, evcent
	#
	# REFERENCES:
	# -> Newman, M.E.J. (2010), "Networks: an introduction", Oxford University Press
	# 
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("igraph")

	if (!is.matrix(sparseA)){
		stop("Input (sparseP) should be a matrix")
	}
	if (!is.matrix(sparseP)){
		stop("Input (sparseP) should be a matrix")
	}
	else if (!isSymmetric(sparseP)){
		stop("Input (sparseP) should be a symmetric matrix")
	}
	else if (!evaluateS(sparseP, verbose = FALSE)$posEigen){
		stop("Input (sparseP) is expected to be positive definite")
	}
	else if (class(as.table) != "logical"){
		stop("Input (as.table) is of wrong class")
	} 
	else{
		# Some warnings
		if (all(sparseA != 0)){warning("Given input (sparseA) implies a saturated conditional independence graph")}
		if (all(sparseA == 0)){warning("Given input (sparseA) implies an empty conditional independence graph")}

		if (all(sparseP != 0)){warning("Given input (sparseP) implies a saturated conditional independence graph")}
		if (all(sparseP[!diag(nrow(sparseP))] == 0)){warning("Given input (sparseP) implies an empty conditional independence graph")}

		###############################################
        	# statistics from A
        	###############################################
	
        	# in and out degree        
        	degreeAout <- ncol(sparseA) - colSums(sparseA==0)
        	degreeAin  <- nrow(sparseA) - rowSums(sparseA==0)
        
		# signs of edges
		nPosAout <- apply(sign(sparseA), 2, function(Z){ sum(Z == 1) }) 
		nNegAout <- apply(sign(sparseA), 2, function(Z){ sum(Z == -1) })

		# signs of edges
		nPosAin <- apply(sign(sparseA), 1, function(Z){ sum(Z == 1) }) 
		nNegAin <- apply(sign(sparseA), 1, function(Z){ sum(Z == -1) })

		# centrality measures of A
	    	# to be included

	    	############################################### 
	    	# statistics from Pe
        	###############################################

        	# (partial) variance of the error
        	pvarsPe <- 1/diag(sparseP)
        	Se      <- solve(sparseP)
        	varsPe  <- diag(Se)

		# signs of edges of P
		slh <- diag(sparseP)
		diag(sparseP) <- 0 
		nPosPe <- apply(sign(sparseP), 2, function(Z){ sum(Z == 1) }) 
		nNegPe <- apply(sign(sparseP), 2, function(Z){ sum(Z == -1) })
        	diag(sparseP) <- slh

	    	# adjacency to graphical object
    		adjMat  <- adjacentMat(sparseP)
		CIGerror <- graph.adjacency(adjMat, mode = "undirected")
	    
		# centrality measures of P
		degreePe          <- degree(CIGerror)
		betweennessPe     <- betweenness(CIGerror)
		closenessPe       <- closeness(CIGerror)
		eigenCentralityPe <- evcent(CIGerror)$vector

		###############################################
		# statistics for Y
		###############################################

		# centrality measures of Py
		CIGy <- graph.adjacency(CIGofVAR1(sparseA, sparseP, "global"), mode="undirected")
		degreePy          <- degree(CIGy)
		betweennessPy     <- betweenness(CIGy)
		closenessPy       <- closeness(CIGy)
		eigenCentralityPy <- evcent(CIGy)$vector
        		
		# calculate the variance of Y	
		Syy <- Se
		for (tau in 1:1000) {
			Atau <- sparseA %^% tau
			Syy <- Syy + Atau %*% Se %*% t(Atau)
			if (max(abs(Atau)) < 10^(-20)){ break }
		}
		varsY <- diag(Syy)

		# Calculate nodes' mutual information
		MI                 <- mutualInfoVAR1(sparseA, sparseP, 2) 
		MIatTplus1         <- MI[,1] 
		names(MIatTplus1)  <- colnames(sparseP)
		MIatTplus2         <- MI[,2] 
		names(MIatTplus2)  <- colnames(sparseP)

		# Calculate nodes' influence response
		IRF                <- impulseResponseVAR1(sparseA, 2)
		IRFatTplus1        <- IRF[,1] 
		names(MIatTplus1)  <- colnames(sparseP)
		IRFatTplus2        <- IRF[,2] 
		names(IRFatTplus2) <- colnames(sparseP)

		# return
		if (as.table){
			networkStats <- cbind(degreeAin, degreeAout, nNegAin, nPosAin, nNegAout, nPosAout, degreePe, betweennessPe, closenessPe, eigenCentralityPe,
                        		nNegPe, nPosPe, varsPe, pvarsPe, varsY, degreePy, betweennessPe, closenessPy, eigenCentralityPy,
                         		MIatTplus1, MIatTplus2, IRFatTplus1, IRFatTplus2)
					colnames(networkStats) <- c("degreeAin", "degreeAout", "nNegAin", "nPosAin", "nNegAout", "nPosAout", "degreePe", 
						"betweennessPe", "closenessPe", "eigenCentralityPe", "nNegPe", "nPosPe", "variancePe", "partialVarPe", "varianceY", 
						"degreePy", "betweennessPy", "closenessPy", "eigenCentralityPy",
                       				"mutualInfo_Tplus1", "mutualInfo_Tplus2", "itemResponse_Tplus1", "itemResponse_Tplus2")
			return(networkStats)
		} 
		if (!as.table){
			return(list(degreeAin=degreeAin, degreeAout=degreeAout, nNegAin=nNegAin, nPosAin=nPosAin, nNegAout=nNegAout, nPosAout=nPosAout, degreePe=degreePe, 
				     betweennessPe=betweennessPe, closenessPe=closenessPe, eigenCentralityPe=eigenCentralityPe, nNegPe=nNegPe, nPosPe=nPosPe, variancePe=varsPe, 
				     partialVarPe=pvarsPe, varianceY=varsY, degreePy=degreePy, betweennessPy=betweennessPy, closenessPy=closenessPy, eigenCentralityPy=eigenCentralityPy,
                     mutualInfo_Tplus1=MIatTplus1, mutualInfo_Tplus2=MIatTplus2, itemResponse_Tplus1=IRFatTplus1, itemResponse_Tplus2=IRFatTplus2))
		}
	}
}




