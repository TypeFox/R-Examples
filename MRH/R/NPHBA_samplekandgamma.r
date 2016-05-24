###########################################################################################
# Created 6Sep13
###########################################################################################

NPHBA_samplekandgamma = function(Mval, Xfixed, delta, outfilename, prune.indc, burnIn, maxIter, thin, Rmp, a,
lambda, betas, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime, checknum, n, dnames,
namesHazGroups, numPHParams, numParams, betaNames, betaLB, betaUB, Xmeans, Xsds, Xstdz, stdIndex,
mat01, inBin, failBin, RmpInd, one_RmpInd, mvec, numHazards, RmpNoPruneIndx, GR, initialValues, RmpRowNums,
indices, continue.chain, k.init, gamma.init, loopctr.start){

    c = d = 1
    if(continue.chain == FALSE){
        gamma.mp = matrix(0.5, nrow = 2^Mval-1, ncol = numHazards)
        k = rep(0.5, numHazards)
        if(GR == TRUE){
            for(rmpctr in 1:numHazards){
                gamma.mp[RmpNoPruneIndx[[rmpctr]],rmpctr] = runif(length(RmpNoPruneIndx[[rmpctr]]), 0, 1)
            }
            k = runif(numHazards, 0, 2)
        }
        write.table(matrix(c('iteration', dnames, paste('beta.', betaNames, sep = ''),
                paste('H00', namesHazGroups, sep = '_'), paste('Rmp', RmpRowNums, sep = ''),
                paste('gammamp', RmpRowNums, sep = ''), paste('k', namesHazGroups, sep = '_'),
                paste('a', namesHazGroups, sep = '_'),paste('lambda', namesHazGroups, sep = '_')), nrow = 1),
                paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
    } else {
        gamma.mp = gamma.init
        k = k.init
    }
    initialValues = c(initialValues, list(gamma.init = gamma.mp, k.init = k))
    # Round values needed: 2^M*numhazard (d) + length(betaNames) + 1*numhazard (H00) +
    #                       (2^M-1)*numhazard (Rmp) + 1*numhazard (a) + 1*numhazard (lambda) +
    #                       (2^M-1)*numhazard (gammamp) + 1*numhazard (k)
    round.values = rep(16, 2^Mval*numHazards + length(betaNames) + numHazards + (2^Mval-1)*numHazards +
                            numHazards + numHazards + (2^Mval-1)*numHazards + numHazards)

	# analyzeIndex holds the columns that need to be analyzed in the burn-in, thinning, and
	# convergence diagnosts/graphs.  Easier than writing out the indices every time those
	# functions are used.
	RmpAnalyzeIndex = 1:((2^Mval-1)*numHazards)
	if(!is.null(prune.indc)){	
		RmpAnalyzeIndex = rep(NA, (2^Mval-1)*numHazards)
		for(rmpctr in 1:numHazards){		
			RmpAnalyzeIndex[1:(2^Mval-1) + (rmpctr-1)*(2^Mval-1)][which(prune.indc[,rmpctr] == 0)] = 
				which(prune.indc[,rmpctr] == 0)+(rmpctr-1)*(2^Mval-1)
		}
	}
	RmpAnalyzeIndex = RmpAnalyzeIndex[!is.na(RmpAnalyzeIndex)]
	analyzeIndex = c(1:numHazards + 2^Mval*(2*numHazards-1) + numPHParams + 1, 
					 RmpAnalyzeIndex + 2^Mval*(2*numHazards-1) + numPHParams + numHazards + 1,
					 RmpAnalyzeIndex + 2^Mval*(2*numHazards-1) + (2^Mval-1)*numHazards + numPHParams + numHazards + 1,
					 1:numHazards + 2^Mval*(2*numHazards-1) + (2^Mval-1)*2*numHazards + numPHParams + numHazards + 1)
	if(numPHParams > 0){	analyzeIndex = c(1:numPHParams + 2^Mval*numHazards + 1, analyzeIndex)	}
	
	#########################################################
	# Run Gibbs sampler
	#########################################################
	outdata = NULL
	loopctr = loopctr.start
	convergence = best.thin.found = FALSE
	F = H = rep(NA, length = n)
	if(fix.thin == TRUE){ best.thin.found = TRUE }
	# Index for convergence graphs and testing
	while((loopctr <= maxIter & convergence == FALSE & fix.max == FALSE) | 
		  (loopctr <= maxIter & fix.max == TRUE)){

		######## Draw H00 ########
		# Calculate F for each hazard group, which is needed in the posterior of each H00
		for(hazCtr in 1:numHazards){
			F[indices[[hazCtr]]] = calc_F(mat01 = mat01, Rmp = Rmp[,hazCtr], M = Mval, 
										  inBinMat = inBin[indices[[hazCtr]],])
		}			
		# Draw H00
		H00 = NULL
		for(hazCtr in 1:numHazards){
			H00 = c(H00, rgamma(1, shape = a[hazCtr]+sum(delta[indices[[hazCtr]]]), 
					scale = (1/lambda[hazCtr]+sum(exp(as.matrix(Xfixed[indices[[hazCtr]],])%*%betas)*
												  F[indices[[hazCtr]]]))^-1))
		}
		######## Draw Rmp values ########
		# Calculate H, which is needed in the posterior of the Rmps
		for(hazCtr in 1:numHazards){
			H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
							inBinMat = inBin[indices[[hazCtr]],])
		}
		# Draw Rmps
		for(hazCtr in 1:numHazards){
			for(rmpctr in RmpNoPruneIndx[[hazCtr]]){
				Rmp[rmpctr, hazCtr]= .Call("arms", c(.001, .999), f = function(x) logRmpPost_nonPHBA(x, 
						RmpFull = Rmp[,hazCtr], H00val = H00[hazCtr], kval = k[hazCtr], aval = a[hazCtr], 
						gamma.mpval = gamma.mp[rmpctr,hazCtr],
						betavec = as.matrix(betas), X.mat = as.matrix(Xfixed[indices[[hazCtr]],]), 
						mval = mvec[rmpctr], RmpIndic = RmpInd[indices[[hazCtr]],rmpctr], 
						one_RmpIndic = one_RmpInd[indices[[hazCtr]],rmpctr],
						deltavals = delta[indices[[hazCtr]]], Mvalue = Mval, 
						inBinMat = inBin[indices[[hazCtr]],], mat01 = mat01, 
						formula = rmpctr), Rmp[rmpctr,hazCtr], as.integer(1), new.env())
			}
		}
		
		##### Draw beta values ########
		# H needs to be recalculated first
		for(hazCtr in 1:numHazards){
			H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00[hazCtr], 
										  Rmp = Rmp[,hazCtr], M = Mval, inBinMat = inBin[indices[[hazCtr]],])
		}
		if(numPHParams > 0){
			if(length(stdIndex) > 0){
				betas.stdz = betas*Xsds
				H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
				for(hazCtr in 1:numHazards){
					H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00.stdz[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
												  inBinMat = inBin[indices[[hazCtr]],])
				}
				for(betaCtr in 1:numPHParams){
					betas.stdz = betas*Xsds
					H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
					for(hazCtr in 1:numHazards){
						H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00.stdz[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
													  inBinMat = inBin[indices[[hazCtr]],])
					}
					betas.stdz[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), f = function(x) 
												logbetaPost_NPHBA(x, betaFull = betas.stdz, deltavals = delta, Xmatrix = Xstdz, Hvals = H, 
																  whichBeta = betaCtr, mu.beta = 0, sigma.beta = 10), betas.stdz[betaCtr], 
												as.integer(1), new.env())
				}
				betas = betas.stdz/Xsds
			} else {
				for(betaCtr in 1:ncol(Xfixed)){
					betas[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), f = function(x) 
										   logbetaPost_NPHBA(x, betaFull = betas, deltavals = delta, Xmatrix = Xfixed, Hvals = H, 
														 whichBeta = betaCtr, mu.beta = 0, sigma.beta = 10), betas[betaCtr], 
										   as.integer(1), new.env())
				}
			}
		}
		
		##### Draw gamma ########
		for(hazCtr in 1:numHazards){
			for(gammactr in RmpNoPruneIndx[[hazCtr]]){
				gamma.mp[gammactr, hazCtr] = .Call("arms", c(.001, .999), 
												   f = function(x) loggammaPost(x, Rmp = Rmp[gammactr,hazCtr], kval = k[hazCtr], aval = a[hazCtr],
																				cval = c, dval = d, mvecval = mvec[gammactr]),
												   gamma.mp[gammactr, hazCtr], as.integer(1), new.env())
			}
		}

		##### Draw k ########
		for(hazCtr in 1:numHazards){
			k[hazCtr] = .Call("arms", c(0, 10), f = function(x) logkPost(x, gamma.mpval = gamma.mp[,hazCtr], a = a[hazCtr], 
																		 Rmpvals = Rmp[,hazCtr], Rmp.exp = 2*gamma.mp[,hazCtr]*a[hazCtr],
																		 one_Rmp.ex = 2*(1-gamma.mp[,hazCtr])*a[hazCtr], 
																		 mvec = mvec, mu.k = 100),
							  k[hazCtr], as.integer(1), new.env())
		}
		
		##### Draw lambda ########
		for(hazCtr in 1:numHazards){
			lambda[hazCtr] = .Call("arms", c(.0001, 10), f = function(x) logLambdaPost(x, 
					mu.lval = 100, H00val = H00[hazCtr], aval = a[hazCtr]), lambda[hazCtr], as.integer(1), new.env())
		}
		
		##### Draw a #########
		for(hazCtr in 1:numHazards){
			a[hazCtr] = sample.a(1:50, lambdaval = lambda[hazCtr], mu.aval = 4, 
								 H00val = H00[hazCtr], kval = k[hazCtr], mvec = mvec, Mval = Mval, 
								 Rmpvec = Rmp[,hazCtr], gamma.mpvec = gamma.mp[,hazCtr])
		}
		
		#######################################################################################
		#					Gather and provide information for user 
		#######################################################################################
        # At 100 iterations, give the user a run time estimate
        if(loopctr == loopctr.start + 99){ checkRunTime(maxIter-loopctr.start+1, systemtime)  }
		# If it is every thin_th iteration, save to the data set
		if((loopctr %% thin) == 0 & loopctr >= burnIn){ 
			d.iter = betaholder = NULL
			for(hazCtr in 1:numHazards){
				d.iter = c(d.iter, calc_dfast(Mval, Rmp[,hazCtr], H00[hazCtr], mat01))
				if(hazCtr > 1){
					betaholder = c(betaholder, log(d.iter[2^Mval*(hazCtr-1)+1:2^Mval]/d.iter[1:2^Mval]))
				}
			}
			if(numPHParams > 0){
				outdata = rbind(outdata, 
								matrix(c(loopctr, d.iter, betas, betaholder, H00,
                                matrix(Rmp, nrow = 1, byrow = TRUE), matrix(gamma.mp, nrow = 1, byrow = TRUE),
                                k, a, lambda), nrow = 1))
			} else {
				outdata = rbind(outdata, 
								matrix(c(loopctr, d.iter, betaholder, H00,
                                matrix(Rmp, nrow = 1, byrow = TRUE), matrix(gamma.mp, nrow = 1, byrow = TRUE),
                                k, a, lambda), nrow = 1))
			}	
		} 
        # If the routine reaches 100,000 iterations, check that the rounding values for the MCMC
        # text file will work (not too large)
        if(loopctr == 100000){
            testData = read.table(paste(outfilename, '/MCMCchains.txt', sep = ''), header = TRUE)
            round.values = FindRoundVals(testData)
        }
        # Every 10,000 or writenum iterations, write out the data set and empty the working one.
        if((loopctr %% writenum) == 0 & loopctr >= burnIn){
            if(nrow(outdata) > 1){
                outdata = cbind(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1])))
            } else {
                outdata = matrix(c(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1]))), nrow = 1)
            }
            write.table(outdata, paste(outfilename, '/MCMCchains.txt', sep = ''),
            col.names = FALSE, row.names = FALSE, append = TRUE)
            if(loopctr != maxIter){	outdata = NULL	}
        }
        # Every 100,000 or checknum iterations, check for convergence with the full data set (minus the burned values)
        if((loopctr %% checknum) == 0 & loopctr >= burnIn){
            chainRes = checkMCMCParams(best.thin.found, analyzeIndex, thin, burnIn, round.values, outfilename,
            convergence, fix.burnIn)
            thin = chainRes$thin
            burnIn = chainRes$burnIn
            convergence = chainRes$convergence
            best.thin.found = chainRes$best.thin.found
        }
		if((loopctr %% 5000) == 0){	print(paste(loopctr, 'iterations completed...'))	}
		
		# Increment loopctr
		loopctr = loopctr+1
	}
	
    # Number of estimated parameters: a, lambda, H, k for each hazard,
    # number of unpruned Rmps and gammas, number of betas
    outputdata = list(assessIndex = analyzeIndex,
                        numberofEstParams = 4*numHazards + 2*length(RmpAnalyzeIndex) + numPHParams,
                        numParams = numParams, failBin = failBin, inBin = inBin, GR = GR, loopctr = loopctr,
                        convergence = convergence, initialValues = initialValues, numPHParams = numPHParams,
                        namesHazGroups = namesHazGroups, burnIn = burnIn, thin = thin, RmpAnalyzeIndex = RmpAnalyzeIndex)
	return(outputdata)
			
}

