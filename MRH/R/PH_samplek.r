###########################################################################################
# This code is used to estimate the parameters of a given data set using the MRH
#	methodology.  One to many covariates can be included under the proportional 
#	hazards assumption.
###########################################################################################

PH_samplek = function(Mval, delta, X, Xstdz, outfilename, burnIn, maxIter, thin, gamma.mp,
GR, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime, checknum, n, betaNames,
numParams, betaLB, betaUB, beta, betas.LB.stdz, betas.UB.stdz, stdIndex, Xmeans, Xsds, mat01,
RmpInd, one_RmpInd, inBin, failBin, mvec, Rmp, a, lambda, RmpNoPruneIndx, RmpRowNums,
initialValues, continue.chain, k.init, loopctr.start){

	#########################################################
	# Start output file for parameter estimates
	#########################################################
	# Write out ds, betas, H00, Rmps
    if(continue.chain == FALSE){
        k = 0.5
        if(GR == TRUE){ k = runif(1, 0, 2)  }
        if(numParams > 0){
            if(length(betaNames) != numParams){	betaNames = 1:numParams	}
                write.table(matrix(c('iteration', paste('d',1:2^Mval, sep = ''), paste('beta', betaNames, sep = '.'),
                    'H00', paste('Rmp', RmpRowNums, sep = ''), 'k', 'a', 'lambda'), nrow = 1),
                    paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
        } else {
                write.table(matrix(c('iteration', paste('d',1:2^Mval, sep = ''),
					'H00', paste('Rmp', RmpRowNums, sep = ''), 'k', 'a', 'lambda'), nrow = 1),
					paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
        }
    } else {
        k = k.init
    }
    initialValues = c(initialValues, list(k.init = k))
    # number of rounding values needed: 2^M (d) + length(betaNames) (optional) + 1 (H00) + 2^M-1 (Rmp) +
    # 1 (a) + 1 (lambda) + k
    round.values = rep(16, 2^Mval + 1 + 2^Mval-1 + 1 + 1 + 1)
    if(numParams > 0){  round.values = c(round.values, rep(16, length(betaNames)))   }
    
	#########################################################
	# Run Gibbs sampler
	#########################################################
	convergence = best.thin.found = FALSE
	if(fix.thin == TRUE){ best.thin.found = TRUE }
	outdata = NULL
	loopctr = loopctr.start
	# assessIndex keeps the index of the parameters that need to be tested for convergence and plotted
	# is equal to the betas, 1 H, 2^M-1 Rmps minus those that are not pruned, 1 k 
	assessIndex = c(1:(numParams+1), (1:(2^Mval-1))[RmpNoPruneIndx]+numParams+1, 2^Mval-1+numParams+1+1)+2^Mval+1

	while((loopctr <= maxIter & convergence == FALSE & fix.max == FALSE) | 
		  (loopctr <= maxIter & fix.max == TRUE)){
		
		######## Draw H00 ########
		# Calculate F, which is needed in the posterior of H00
		F = calc_F(mat01 = mat01, Rmp = Rmp, M = Mval, inBinMat = inBin)
		# Draw H00
		H00 = rgamma(1, shape = a+sum(delta), scale = (1/lambda+sum(exp(X%*%beta)*F))^-1)
	
		######## Draw Rmp values ########
		# Calculate H, which is needed in the posterior of the Rmps
		H = calc_H(mat01 = mat01, H00 = H00, Rmp = Rmp, M = Mval, inBinMat = inBin)
		# Draw Rmps
		for(Rmpctr in RmpNoPruneIndx){
			Rmp[Rmpctr] = .Call("arms", c(.001, .999), f = function(x) logRmpPost_PHcovs(x, RmpFull = Rmp, H00val = H00, kval = k, 
																						 aval = a, 
																						 gamma.mpval = gamma.mp[Rmpctr],
																						 betavec = beta, X.mat = X, 
																						 mval = mvec[Rmpctr], 
																						 RmpIndic = RmpInd[,Rmpctr], 
																						 one_RmpIndic = one_RmpInd[,Rmpctr],
																						 deltavals = delta, Mvalue = Mval, 
																						 inBinMat = inBin, mat01 = mat01, 
																						 formula = Rmpctr), Rmp[Rmpctr], 
								as.integer(1), new.env())
		}
		
		##### Draw beta ########
		# H needs to be recalculated first
		if(numParams > 0){
			if(length(stdIndex) > 0){
				betas.stdz = beta*Xsds
				H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
				H = calc_H(mat01 = mat01, H00 = H00.stdz, Rmp = Rmp, M = Mval, inBinMat = inBin)
				for(betaCtr in 1:numParams){
					H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
					H = calc_H(mat01 = mat01, H00 = H00.stdz, Rmp = Rmp, M = Mval, inBinMat = inBin)
					betas.stdz[betaCtr] = .Call("arms", c(betas.LB.stdz[betaCtr], betas.UB.stdz[betaCtr]), 
													   f = function(x) logbetaPost_PH(x, betaFull = betas.stdz, 
																					  deltavals = delta, Xmatrix = Xstdz, 
																					  Hvals = H, whichBeta = betaCtr, mu.beta = 0, 
																					  sigma.beta = 10), betas.stdz[betaCtr], 
													   as.integer(1), new.env())
				}
				beta = betas.stdz/Xsds
			} else {
				H = calc_H(mat01 = mat01, H00 = H00, Rmp = Rmp, M = Mval, inBinMat = inBin)
				for(betaCtr in 1:numParams){
					beta[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), 
										  f = function(x) logbetaPost_PH(x, betaFull = beta, deltavals = delta, 
																		 Xmatrix = X, Hvals = H, whichBeta = betaCtr, mu.beta = 0, 
																		 sigma.beta = 10), beta[betaCtr], 
										  as.integer(1), new.env())
				}
			}
		}
		##### Draw k ########
		k = .Call("arms", c(0, 10), f = function(x) logkPost(x, gamma.mpval = gamma.mp, a = a, Rmpvals = Rmp,
																Rmp.exp = 2*gamma.mp*a,
																one_Rmp.ex = 2*(1-gamma.mp)*a, mu.k = 100, mvec = mvec),
				  k, as.integer(1), new.env())
		
		##### Draw lambda ########
		lambda = .Call("arms", c(.0001, 10), f = function(x) logLambdaPost(x, mu.lval = 100, H00val = H00, aval = a),
					   lambda, as.integer(1), new.env())
		
		##### Draw a #########
		a = sample.a(1:50, lambdaval = lambda, mu.aval = 4, H00val = H00, kval = k, 
							 mvec = mvec, Mval = Mval, Rmpvec = Rmp, gamma.mpvec = gamma.mp)

        #######################################################################################
        #					Gather and provide information for user
        #######################################################################################
        # At 100 iterations, give the user a run time estimate
        if(loopctr == loopctr.start + 99){ checkRunTime(maxIter-loopctr.start+1, systemtime)  }
        # Every thinth iteration after the burn-in number, save the data set.
        if((loopctr %% thin) == 0 & loopctr >= burnIn){
            d.iter = calc_dfast(Mval, Rmp, H00, mat01)
            if(numParams > 0){
                outdata = rbind(outdata, matrix(c(loopctr, d.iter, beta, H00, Rmp, k, a, lambda), nrow = 1))
            } else {
                outdata = rbind(outdata, matrix(c(loopctr, d.iter, H00, Rmp, k, a, lambda), nrow = 1))
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
                outdata = matrix(c(outdata[,1],
                            sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1]))), nrow = 1)
            }
            write.table(outdata, paste(outfilename, '/MCMCchains.txt', sep = ''),
                        col.names = FALSE, row.names = FALSE, append = TRUE)
            if(loopctr != maxIter){	outdata = NULL	}
        }
        # Every 100,000 or checknum iterations, check for convergence with the full data set (minus the burned values)
        if((loopctr %% checknum) == 0 & loopctr >= burnIn){
            chainRes = checkMCMCParams(best.thin.found, assessIndex, thin, burnIn, round.values, outfilename,
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

	outputdata = list(assessIndex = assessIndex, numberofEstParams = 2 + 1 + numParams + length(RmpNoPruneIndx) + 1,
                        numParams = numParams, failBin = failBin, inBin = inBin, GR = GR, loopctr = loopctr,
                        convergence = convergence, initialValues = initialValues, burnIn = burnIn, thin = thin)
	return(outputdata)
	
}
