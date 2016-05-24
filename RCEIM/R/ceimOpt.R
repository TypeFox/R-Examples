ceimOpt <- function(OptimFunction="testFunOptimization", nParams=1, 
				   minimize=TRUE,
				   Ntot=1000, N_elite=floor(Ntot/4), N_super=1, alpha=1, epsilon=0.1, q=2,
				   maxIter=50, waitGen=maxIter, boundaries=t(matrix(rep(c(-10,10),nParams), ncol=nParams)), plotConvergence=FALSE, 
				   chaosGen=maxIter,
				   handIterative=FALSE, verbose=FALSE, plotResultDistribution=FALSE,
				   parallelVersion=FALSE) {
#
# This is the Cross-Entropy Inspired Method implementation for optimization
#
# "OptimFunction" is a string with the name of the function that will be optimized
#
# "minimize" is a boolean indicating if the OptimFunction will be minimized or maximized
#
# "nParams" is an integer with the number of parameters
#
# "Ntot" is an integer with the number of individuals per iteration
#
# "N_elite" is an integer with the number of elite individuals, or in other words, the individuals 
#           used to define the individuals of the new iteration
#
# "N_super" is an integer with the number of super-individuals, or those individuals with the 
#           best fitness value, that are copied to the new iteration
# 
# "alpha, q" are parameters of the CE method used to control the convergence rate, and to prevent
#            early convergence to local optima
#
# "epsilon" is a convergence control parameter: if the maximum st.dev. of the parameters of the 
#           elite individuals divided by its average value is smaller than this number, the method 
#           considers that it converged
#
# "maxIter" is the absolute maximum number of iterations that the method will run
#
# "waitGen" is the number of iterations that the method will wait: after "waitGen" without any
#           improvement in the best individual, the method gives up and return the best individual
#           as an answer
#
# "Boundaries" is a matrix with as many rows as there are parameters and two columns
#            the first column stores the minimum value, while the second, the maximum
#
# "chaosGen" is the number of iterations before the method replaces all the solutions, but the 
#            super-individuals, by a new random trial
#
# "plotConvergence" is a flag to indicate if the user wants to check visually the convergence of 
#            the method
#
# "handIterative" is a flag to indicate if the user wants to press enter between the iterations
# 
# "verbose" is a flag to indicate if the user wants to receive some convergence and distribution
#            information printed on the screen
# 
# "parallelVersion" is a flag to indicate if the user wants to use all the cores in his/her computer
#            to compute the fitnees functions
#
#				   maxIter=50, waitGen=maxIter, boundaries=t(as.matrix(c(-10, 10))), plotConvergence=FALSE,
	
	# ---- Initialization --------------------------------------------------------------------------
	# Check if the boundaries are compatible with nParams
	if(dim(boundaries)[1] != nParams) {
		stop("Boundaries dimensions (n. rows) are not compatible with the number of parameters.")
	}	
	
	# Check if boundaries are valid
	if(dim(boundaries)[2] != 2) {
		stop("Boundaries dimensions (n. cols) are not valid. It should be two (min, max)!")
	}
	if(sum(boundaries[,1] < boundaries[,2]) != length(boundaries[,1])) {
		stop("Trying to use invalid boundaries, with Min > Max!")
	}
	
	# Should we plot the convergence
	if(plotConvergence) { dev.new() }
	
	# Should we minimize or maximize?
	if(minimize) { mfactor <- 1 } 
	else { mfactor <- -1 }
	# ----------------------------------------------------------------------------------------------
		
	# Declare some vectors...
	convergenceStatsBest <- vector("double", maxIter)
	convergenceStatsMean <- vector("double", maxIter)
	convergenceStatsSdev <- vector("double", maxIter)
	
	Sfunc <- match.fun(OptimFunction)
	
	# Generate a uniformly random first generation within the boudaries
	# After all, at this first iteration we can pretend to have no information whatsoever
	# about our problem. 
	param_Estimated <- matrix(0, nrow=Ntot, ncol=nParams)
	for(i in 1:nParams) {
		param_Estimated[,i] <- runif(Ntot, min=boundaries[i,1], max=boundaries[i,2])
	}
	
	# Compute the initial statistics
	mu <- vector("double", nParams)
	sig <- vector("double", nParams)
	for(i in 1:nParams) {
		mu[i] <- mean(param_Estimated[,i])
		sig[i] <- sd(param_Estimated[,i])
	}
	fracsig <- 10 * epsilon

	# Declare some variables	
	iter <- 1
	waitGenCounter <- 0
	chaosCounter <- 0
	oldEliteS <- 0
	Svals <- vector("double", Ntot)
	Crit <- "Temp"
	
	# now run the internal main loop, 
	# while the iteration number is smaller than maxIter
	# or stop due to some other criteria, such as the fractional epsilon value
	#while( (iter < maxIter) && (epsilon < max(sig) ) && (waitGenCounter <= waitGen) && (Crit!="SAME") ) {	
	while( (iter < maxIter) && (epsilon < fracsig ) && (waitGenCounter <= waitGen) && (Crit!="SAME") ) {	
		
		# Print some convergence information...
		if(verbose) {
			for(i in 1:nParams) {
				#cat(paste("[",iter," - ", i,"] - ", round(mu[i],3) ," - ", round(sig[i],3), "\n"))
				cat(paste("[",iter," - ", i,"] - ", round(mu[i],4) ," - ", round(fracsig,4), "\n"))
			}
		}
		
		# Now compute S (the fitness) for each individual
		if(parallelVersion){
      if(requireNamespace("parallel", quietly=TRUE)) {
			  #library(parallel)
			  pp <- parallel::mclapply(1:Ntot, function(x) {Sfunc(param_Estimated[x,])})						         
			  for(i in 1:Ntot) {
				  #Svals[i] <- mfactor * Sfunc(param_Estimated[i,])
			  	Svals[i] <- mfactor * pp[[i]]
			  }
      }
		} else {
			for(i in 1:Ntot) {
				Svals[i] <- mfactor * Sfunc(param_Estimated[i,])
			}
		}
		
		# Get the N_elite individuals with the best S values
		dfsel <- data.frame(param_Estimated, S=Svals)
		dfsel <- sortDataFrame(dfsel, "S")
		sN <- length(dfsel$S)
		elite <- dfsel[1:N_elite,]
		if(verbose) {
			cat(paste("[",iter,"      ] -                       ", round(elite[1,(length(elite))],3), "\n"))
		}
		
		convergenceStatsBest[iter] <- elite$S[1]
		convergenceStatsMean[iter] <- mean(elite$S)
		convergenceStatsSdev[iter] <- sd(elite$S)
				
		# Check if the best individual changed or not
		if (abs(oldEliteS - elite$S[1]) > 0) {
			# if it did, put a zero at the waitGenCounter
			waitGenCounter <- 0
			chaosGenCounter <- 0
		} else {
			# if not, add one to the waitGenCounter
			waitGenCounter <- waitGenCounter + 1
			chaosGenCounter <- chaosGenCounter + 1
		}
		oldEliteS <- elite$S[1]

		# Check if the mean and the best individual have the same value		
		if(verbose) {
				cat(paste("[",iter," - ", maxIter,"] -  CONVERGENCE STATUS \t\t\t ", round(convergenceStatsBest[iter],3) ," - ", round(convergenceStatsMean[iter],3)," - ", round(convergenceStatsSdev[iter],3), "\n"))
		}
		if(convergenceStatsBest[iter] == convergenceStatsMean[iter]) {
			Crit <- "SAME"
		}
		
		# Now, generate the next generation of individuals...
		if(epsilon < min(sig) ) {
			alpha_d <- alpha - alpha*(1-1/iter)^q
			for(i in 1:nParams) {
				mu[i] <- alpha * mean(elite[,i]) + (1-alpha) * mu[i] 
				sig[i] <- alpha_d * sd(elite[,i]) + (1-alpha_d) * sig[i]
			}
			param_Estimated <- matrix(0, nrow=Ntot, ncol=nParams)
			Nnew <- Ntot-N_super
			for(i in 1:nParams) {
				param_Estimated[1:Nnew,i] <- rnorm(Nnew, mu[i], sig[i])
			}
			# If the user wants to keep super-individuals, add them to the list now...
			if(N_super > 0) {
				for(i in 1:nParams) {
					param_Estimated[Ntot:(Nnew+1),i] <- elite[1:N_super,i]
				}
			}
			maxSigIdx <- which(sig==max(sig))
			fracsig <- sig[maxSigIdx]/abs(mu[maxSigIdx])
		}
		
		# if the user wants to add chaos after a certain number of generations without improvement...
		if(chaosGenCounter >= chaosGen) {
			if(verbose) {
				cat(paste("[",iter," - ", i,"] -  Adding CHAOS! \n"))
			}
			for(i in 1:nParams) {
#				param_Estimated[1:Nnew,i] <- runif(Nnew, min=boundaries[i,1], max=boundaries[i,2])               # uniformly random
#				param_Estimated[1:Nnew,i] <- rnorm(Nnew, mean=elite[i,1], sd=(boundaries[i,2]-boundaries[i,1]) ) # normaly distributed, taking the best as the mean
				param_Estimated[1:Nnew,i] <- rnorm(Nnew, mean=mu[i], sd=(boundaries[i,2]-boundaries[i,1])/2 )    # normaly distributed, taking the average as the mean
			}
			if(N_super > 0) {
				for(i in 1:nParams) {
					param_Estimated[Ntot:(Nnew+1),i] <- elite[1:N_super,i]
				}
			}
			chaosGenCounter <- 0
		}
		param_Estimated <- enforceDomainOnParameters(param_Estimated, boundaries)
		
		# if the user wants to check visually the convergence of the problem
		if(plotConvergence) {
			plot(1:iter, convergenceStatsMean[1:iter], type="n", xlab="Iteration", ylab="Fval", 
				ylim=c(min(c(convergenceStatsBest, convergenceStatsMean)),max(convergenceStatsMean)), main="Convergence status")
				#ylim=c(min(c(convergenceStatsBest*0.9, convergenceStatsMean*0.9)),max(convergenceStatsMean)), main="Convergence status")
			    # , log="y")
			grid()
			#overPlotErrorPolygon(1:iter, convergenceStatsMean[1:iter],   convergenceStatsSdev[1:iter], col=rgb(0,0,1,0.25), logPlot=TRUE, border=NA)
			overPlotErrorPolygon(1:iter, convergenceStatsMean[1:iter],   convergenceStatsSdev[1:iter], col=rgb(0,0,1,0.25), logPlot=FALSE, border=NA)
			lines(1:iter, convergenceStatsMean[1:iter], col="black", lwd=2)
			lines(1:iter, convergenceStatsBest[1:iter], col="red", lwd=2)			
		}
		if(handIterative) {
			cat(" Press enter for the next iteration... ")						
			readline()
		}
		iter <- iter + 1
	}
	
	# And finally, the convergence informations...
	if(iter >= maxIter) { 
		Convergence <- FALSE
		Crit <- "Maximum iterations reached."
	} else if (waitGenCounter >= waitGen) {
		Convergence <- FALSE
		Crit <- "Maximum wait generations reached without changes in the best individual."
	} else if(Crit=="SAME") {
		Convergence <- TRUE
		Crit <- "Convergence :: (Mean-Best) are the same, and the stdev is small."
	} else {
		Convergence <- TRUE
		Crit <- "Convergence criteria reached :: Max(sigma/mu) is smaller than the requested value."
	}

	
	if(plotResultDistribution) {
		plotEliteDistrib(elite)
	}
	
	# We are returning the best individual. Some people believe that the average should be returned... 
	# why? simple, with the average you have a estimate of the error (the st.dev.), while with the best
	# you don't... but you know that the best individual is "better" than the average, and we want 
	# the best solution we can provide. Now one should think about how to estimate an error for this 
	# "best" individual.
	return(list(BestMember=unlist(elite[1,1:(nParams+1)]), Convergence=Convergence, Criteria=Crit, Iterations=iter, EliteMembers=elite))
}