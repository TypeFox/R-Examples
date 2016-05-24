#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: phase3.r
# *
# * Description: This module contains the function phase3 which runs final
# * iterations
# * with constant theta for estimating rate parameter, derivative matrix and
# * covariance matrix
# *****************************************************************************/
##args: x model object, z control object
##@phase3 siena07 Does phase 3
phase3 <- function(z, x, ...)
{
    ## initialize phase 3
    f <- FRANstore()
    DisplayTheta(z)
    z$Phase <-  3
    int <- z$int
    z$returnDeps <- z$returnDepsStored

    if (x$checktime) z$ctime <- proc.time()[3]

    ## fix up iteration numbers if using multiple processors
    if (x$n3 %% int > 0)
    {
        endNit <- x$n3 + int - x$n3 %%int
    }
    else
    {
        endNit <- x$n3
    }
	z$n3 <- endNit
    z$Phase3nits <- endNit

    ## revert to original requested method for phase 3 unless symmetric
    if (z$FinDiff.method && !x$FinDiff.method &&
        (!is.null(z$FinDiffBecauseSymmetric)) && z$FinDiffBecauseSymmetric)
    {
        z$Deriv <- FALSE
        z$FinDiff.method <- TRUE
    }
    else
    {
        z$Deriv <- !x$FinDiff.method
        z$FinDiff.method <- x$FinDiff.method
    }
    if (z$FinDiff.method)
        Report('Estimation of derivatives by the finite difference method.\n\n',
               outf)
    else
        Report('Estimation of derivatives by the LR method (type 1).\n\n', outf)

	z <- createSiena07stores(z, z$n3, f)
    xsmall<- NULL
    zsmall <- makeZsmall(z)
    if (!x$maxlike && !is.null(z$writefreq))
    {
        if (z$FinDiff.method)
            z$writefreq <- z$writefreq %/% z$pp
        else
            z$writefreq <- z$writefreq %/% 2
        z$writefreq <- roundfreq(z$writefreq)
    }
    z <- AnnouncePhase(z, x)
    Report('Simulated values, phase 3.\n', cf)
    nits <- seq(1, endNit, int)
    nits11 <- min(c(endNit, nits[nits >= 11]))
    writefreq <- z$writefreq
    if (is.null(z$writefreq))
    {
        z$writefreq <- 10
    }
	z <- doPhase1or3Iterations(3, z, x, zsmall, xsmall, nits, 0, nits11,
							   writefreq)
    z <- phase3.2(z,x)
    z
}

##@phase3.2 siena07 Processing at end of phase 3
phase3.2 <- function(z, x, ...)
{
    z$timePhase3 <- (proc.time()['elapsed'] - z$ctime) / z$Phase3nits
    if (x$checktime)
	{
        Report(c('Time per iteration in phase 3   = ',
                 format(z$timePhase3, nsmall=4, digits=4), '\n'), lf)
	}
	z <- CalculateDerivative3(z, x)
	if (!x$simOnly)
	{
		z <- PotentialNR(z, x, FALSE)
	}
    if (any(z$newfixed))
    {
        Report('There was a problem in obtaining convergence.\n', outf)
        Report(c('Therefore, the program decided tentatively to fix',
				 'parameter(s)',
               c(1:z$pp)[z$newfixed], '.\n'), outf)
        Report(c('It may be better to start all over again,',
                 'with better initial values or a reduced model.\n',
                 '(Check that you entered the data properly!)\n'), outf)
    }
    Heading(2, outf, c('End of stochastic approximation algorithm, phase ',
                       z$Phase, '.'))
    Report(c('Total of', z$n,'iterations.\n'), outf)
    Report(c('Parameter estimates based on', z$n - z$Phase3nits,
             'iterations,\n'), outf)
    if (!x$maxlike && z$cconditional)
	{
        Report(c('basic rate parameter',
                 c('', 's')[as.integer(z$f$observations > 2) + 1],
                 ' as well as \n'), sep='', outf)
	}
	Report(c('convergence diagnostics, covariance and derivative matrices ',
			 'based on ', z$Phase3nits, ' iterations.\n\n'), sep='', outf)
    Report('Information for convergence diagnosis.\n', outf)
    Report(c('Averages, standard deviations, ',
           'and t-ratios for deviations from targets:\n'), sep='', outf)
  #  Report(c(date(),'\n'),bof)
    if (x$maxlike)
	{
        Report('\nMaximum Likelihood estimation.', bof)
	}
    else
	{
		if (z$cconditional)
		{
			Report('\nconditional moment estimation.', bof)

		}
		else
		{
			Report('\nunconditional moment estimation.', bof)
		}
	}
    Report('\nInformation for convergence diagnosis.\n', bof)
    Report(c('Averages, standard deviations, ',
        'and t-ratios for deviations from targets:\n'), bof, sep='')
    ##calculate t-ratios
    dmsf <- diag(z$msf)
    sf <- colMeans(z$sf)
	# TS: I wonder why "use" and "use2" are the names in the following lines;
	# the coordinates with "use" are <<not>> used.
    use <- dmsf < 1e-20 * z$scale * z$scale
    use2 <- abs(sf) < 1e-10 * z$scale
    dmsf[use] <- 1e-20 * z$scale[use] * z$scale[use]
    tstat <- rep(NA, z$pp)
    tstat[!use]<- sf[!use] / sqrt(dmsf[!use])
    tstat[use & use2] <- 0
    tstat[use & !use2] <- 999
    z$tstat <- tstat
	# tconv.max = Maximum value of t-ratio for convergence,
	# for any linear combination.
	z$tconv.max <- NA
	if (sum(!z$fixed) > 0)
	{
		mean.dev <- colSums(z$sf)[!z$fixed]/dim(z$sf)[1]
		cov.dev <- z$msf[!z$fixed,!z$fixed]
		if (inherits(try(thisproduct <- solve(cov.dev, mean.dev)),"try-error"))
		{
			Report('Maximum t-ratio for convergence not computable.\n', outf)
		}
		else
		{
			z$tconv.max <- sqrt(t(mean.dev) %*% thisproduct)
		}
	}
    mymess1 <- paste(format(1:z$pp,width=3), '. ',
                    format(round(sf, 4), width=8, nsmall=4), ' ',
                    format(round(sqrt(dmsf), 4) ,width=8, nsmall=4), ' ',
                    format(round(tstat, 4), width=8, nsmall=3), sep='')
    mymess2 <- c('', '    (fixed parameter)')[as.numeric(z$fixed) + 1]
    mymess <- paste(mymess1, mymess2)
    PrtOutMat(as.matrix(mymess), outf)
    PrtOutMat(as.matrix(mymess1), bof)
    ##  Report(mymess1, bof, fill=80)
    tmax <- max(abs(tstat)[!z$fixed & !z$BasicRateFunction & z$resist > 0.9])
    z$tconv <- tstat
    error <- (abs(tmax) > 4.0 / sqrt(z$Phase3nits)) && (abs(tmax) > 0.3)
    if (tmax >= 0.4 & !z$error)
	{
        z$error <- TRUE
	}
	Report('Good convergence is indicated by the t-ratios ', outf)
    if (any(z$fixed))
	{
		Report('of non-fixed parameters ', outf)
	}
    Report('being close to zero.\n', outf)
    if (z$Phase3nits < 100)
	{
        Report(c('(Since the diagnostic checks now are based only on ',
                 z$Phase3nits,
                 ' iterations,', '\nthey are not reliable.)\n'), sep='', outf)
	}
    if (error) ## also test subphase here but not relevant to phase 3, I think
    {
        Report('One or more of the t-statistics are rather large.\n', outf)
        if (tmax > 0.5)
		{
            Report('Convergence of the algorithm is doubtful.\n', outf)
		}
		## removed repfortotal loop possibility here as not functioning now
        if (z$Phase3nits <= 50)
		{
            Report(c('Note that the standard deviations are based on',
                     'few simulations.\n'), outf)
		}
	}
    if (x$maxlike)
    {
        Report('Autocorrelations during phase 3 : \n', outf)
        Report(paste(format(1:z$pp, width=3), '. ',
                     format(z$sfl, width=8, digits=4),
                     '\n', collapse="", sep=""), outf)
        Report ('\n', outf)
        Report('Autocorrelations during phase 3 : \n', cf)
        Report(paste(format(1:z$pp, width=3), '. ',
                     format(z$sfl, width=8, digits=4),
                     '\n', collapse="", sep=""), cf)
        Report ('\n', cf)
    }
    for (j in 1:z$pp)
	{
        if (z$diver[j]) ### don't understand this condition, as AllFixed is true
        {
            Report(c('Warning. Extremely large standard error of parameter',j,
                     '.\n'), outf)
            if (sf[j] < 0.5 * sqrt(dmsf[j]))
			{
                Report('Presumably this parameter must be fixed.\n', outf)
			}
			else
			{
                Report('Maybe the algorithm diverged.\n', outf)
			}
		}
	}
	if (!x$simOnly)
	{
		if (x$maxlike)
		{
			Report('Estimated complete data information matrix: \n', cf)
			PrtOutMat(z$dfra, cf)
			Report(c('Estimated conditional covariance matrix score function ',
				'(unobserved information):\n'), cf)
			PrtOutMat(z$msf, cf)
			Report('\n', cf)
			dfrac <- z$dfra - z$msf
			## dfrac[z$fixed[row(dfrac)] | z$fixed[col(dfrac)]] <- 0
			## a clever way to do it
			dfrac[z$fixed, ] <- 0
			dfrac[ ,z$fixed] <- 0
			diag(dfrac)[z$fixed] <- 1
			if (inherits(try(cov <- solve(dfrac)),"try-error"))
			{
				Report('Noninvertible estimated covariance matrix : \n', outf)
				cov <- NULL
			}
		}
		else
		{
			cov <- z$dinv %*% z$msfc %*% t(z$dinv)
		}
		error <- FALSE
		if (inherits(try(msfinv <- solve(z$msfc)), "try-error"))
		{
			Report('Covariance matrix not positive definite: \n', outf)
			if (any(z$fixed || any(z$newfixed)))
			{
				Report(c('(This may be unimportant, and related to the fact\n',
					'that some parameters are fixed.)\n'), outf)
			}
			else
			{
				Report(c('This may mean that the reported standard errors ',
						'are invalid.\n'), outf)
			}
			z$msfinv <- NULL
		}
		else
		{
			z$msfinv <- msfinv
		}
		if (!is.null(cov))
		{
			z$diver <- (z$fixed | z$diver | diag(cov) < 1e-9) & (!z$AllUserFixed)
			## beware: recycling works for one direction but not the other
			diag(cov)[z$diver] <- 99 * 99
			cov[z$diver, ] <- rep(Root(diag(cov)), each=sum(z$diver)) * 33
			diag(cov)[z$diver] <- 99 * 99
			cov[, z$diver] <- rep(Root(diag(cov)), sum(z$diver)) * 33
			diag(cov)[z$diver] <- 99 * 99
		}
		z$covtheta <- cov
	}
	## ans<-InstabilityAnalysis(z)
	z
}

##@CalculateDerivative3 siena07 Calculates derivative at end of phase 3
CalculateDerivative3<- function(z,x)
{
    z$mnfra <- colMeans(z$sf)
    if (z$FinDiff.method || x$maxlike)
    {
		dfra <- t(as.matrix(Reduce("+", z$sdf) / length(z$sdf)))
		z$regrCoef <- rep(0, z$pp)
		z$regrCor <- rep(0, z$pp)
    }
	else
    {
        dfra <-  derivativeFromScoresAndDeviations(z$ssc, z$sf2)
        if (any(diag(dfra) < 0))
        {
            sub <- which(diag(dfra) < 0)
            dfra[sub,] <- 0
            dfra[,sub] <- 0
            dfra[sub, sub] <- 1
            Report(c("Warning: diagonal element(s)", sub,
                     " of derivative matrix < 0\n"), cf)
        }
		scores <- apply(z$ssc, c(1,3), mean)
		for (i in 1:z$pp)
		{
			if (var(scores[,i]) > 0)
			{
				z$regrCoef[i] <- cov(z$sf[,i], scores[,i])/var(scores[,i])
				z$regrCor[i] <- cor(z$sf[,i], scores[,i])
			}
		}
		Report('Correlations between scores and statistics:\n', cf)
		PrtOutMat(format(as.matrix(t(z$regrCor)), digits = 2, nsmall = 2), cf)
    }
    z$diver <- rep(FALSE, z$pp)
    if (z$AllUserFixed & any(abs(diag(dfra)) < 1e-6))
	{
        z$diver[abs(diag(dfra)) < 1e-6] <- TRUE
	}
	z$msf <- cov(z$sf)
    if (z$Phase3nits > 2)
    {
        z$sfl <- apply(z$sf, 2,
                       function(x)acf(x, plot=FALSE, lag.max=1)[[1]][[2]])
    }
    z$dfra1 <- z$dfra
    z$dfra <- dfra
    z
}

##@PotentialNR siena07 Calculates change if NR step done now
PotentialNR <-function(z,x,MakeStep=FALSE)
{
    z$dfrac <- z$dfra
    z$msfc <- z$msf
    if (!z$AllUserFixed)
    {
        z$dfrac[z$fixed, ] <- 0
        z$dfrac[, z$fixed] <- 0
        diag(z$dfrac)[z$fixed]<- 1
        z$msfc[z$fixed, ] <- 0
        z$msfc[, z$fixed] <- 0
        diag(z$msfc)[z$fixed] <- 1
    }
    if (inherits(try(dinv <- solve(z$dfrac)), "try-error"))
    {
        Report('Error message from inversion of dfra: \n', cf)
        diag(z$dfrac) <- diag(z$dfrac)+0.1*z$scale
        Report('Intervention 3.4: ridge added after phase 3.\n', cf)
        if (inherits(try(dinv <- solve(z$dfrac)), "try-error"))
        {
            Report(c('Warning. After phase 3, derivative matrix non-invertible',
                     'even with a ridge.\n'), cf)
            fchange <- 0
            z$dinv <- z$dfrac
			z$dinv[,] <- 999
        }
        else
        {
            fchange <- dinv %*% colMeans(z$sf)
            z$dinv <- dinv
        }
    }
    else
    {
        fchange <- dinv%*%colMeans(z$sf)
        z$dinv <- dinv
    }
	z$fchange <- fchange
    Report('dfrac :\n', cf)
    PrtOutMat(z$dfrac, cf)
    Report('inverse of dfra :\n', cf)
    PrtOutMat(z$dinv, cf)
    Report(c('A full Quasi-Newton-Raphson step after phase 3\n',
             'would add the following numbers to the parameters, yielding ',
             'the following results:\n'), sep='', cf)
    Report('         change     new value \n', cf)
    Report(c(paste('  ', format(1:z$pp, width=2), '. ',
                   format(round(-fchange, digits=6), width=12, nsmall=6),
                   format(round(z$theta - fchange, 6), width=12, nsmall=6),
                   sep='', collapse='\n'), '\n'), cf)
    if (MakeStep) ##currently not used
    {
        Report(c('\nAt the end of phase ', z$phase,
				 ', parameter values are \n'), outf)
        Report(paste(1:z$pp,'. ',format(z$theta,width=18,digits=6)),outf)
        Report(c('A full Quasi-Newton-Raphson step after phase 3\n',
				 'would add the ',
                 'following numbers to the parameters:\n'), outf)
        Report(paste(1:z$pp, '. ', format(-round(fchange, 6), width=12)), outf)
        Report('\n\n', outf)
        if (z$SomeFixed)
		{
			fchange[z$fixed] <- 0
		}
        ##check positivity
        if (z$anyposj)
        {
            neg <- z$posj & fchange >= z$theta
            fchange[neg] <- 0.5 * z$theta[neg]
            Report(c('Intervention 3.4.3: positivity restriction after phase 3',
                     'coordinate(s)', neg, '.\n'), cf)
        }
    }
    z
}

doPhase1or3Iterations <- function(phase, z, x, zsmall, xsmall, nits, nits6=0,
								  nits11=0, writefreq)
{
	int <- z$int
	for (nit in nits)
	{
		z$nit <- nit
		## initial steps differ between phases
		if (phase == 1)
		{
			if (length(nits) > 1 && nit == nits[2])
			{
				time1 <- proc.time()['elapsed']
				if (x$checktime)
				{
					z$ctime <- time1
				}
			}
			if (nit == nits6)
			{
				time1 <- proc.time()['elapsed']-time1
				if (time1 > 1e-5)
				{
					z$writefreq <- round(20.0/time1)
				}
				else
				{
					z$writefreq <- 5
				}
				if (is.batch())
				{
					z$writefreq <- 10 * z$writefreq
				}
				z$writefreq <- roundfreq(z$writefreq)
			}
			if (nit > nits6 && is.null(z$ctime))
			{
				time1 <- proc.time()['elapsed']
				if (x$checktime)
				{
					z$ctime <- time1
				}
			}
			DisplayIteration(z)
		}
		else ## phase 3
		{
			if (is.null(writefreq))
			{
				if (nit == nits[2])
				{
					z$time1 <- proc.time()[[3]]
				}
				else if (nit == nits11)
				{
					time1 <- proc.time()[[3]] - z$time1
					if (time1 > 1e-5)
					{
						z$writefreq <- max(1, round(20.0 / time1))
					}
					else
					{
						z$writefreq <- 20
					}
					if (is.batch())
					{
						z$writefreq <-  z$writefreq * 2 ##compensation for it
						## running faster with no tcl/tk
					}
					z$writefreq <- roundfreq(z$writefreq)
					writefreq <- z$writefreq
				}
			}
			if (nit <= 5 || nit == 10 || (int==1 && nit %% z$writefreq == 0 ) ||
				(int > 1 && (z$writefreq + 1) < z$n3%/%int &&
				 nit %in% nits[seq(z$writefreq + 1, x$n3 %/% int,
								   z$writefreq)]))
			{
				if (!is.batch())
				{
					DisplayIteration(z)
					##  DisplayTheta(z)
					DisplayDeviations(z, z$sf[nit-1,])
				}
				## else
				## {
				##     Report(c('Phase ', z$Phase,' Iteration ',nit,'\n'))
				## }
				if (nit %% z$writefreq == 0 || (int > 1 &&
						   nit %% z$writefreq == 1))
				{
					increment <- ifelse(nit <= 5, int,
										ifelse(nit <= 10, 5, z$writefreq * int))
					val<- getProgressBar(z$pb)
					if (z$FinDiff.method)
						val <- val + increment * (z$pp + 1)
					else
						val <- val + increment
					## sort out progress bar
					z$pb <- setProgressBar(z$pb, val)
					if (is.batch())
						Report(c("Phase ", z$Phase, " Iteration ", nit,
								 " Progress ",
								 round(val / z$pb$pbmax * 100), "%\n"), sep='')
				}
			}
		}
		## iteration proper
		if (z$int == 1)
		{
			zz <- x$FRAN(zsmall, xsmall)
			if (!zz$OK)
			{
				z$OK <- zz$OK
				z$zz <- zz
				return(z)
			}
			z$n <- z$n + 1
			if (phase == 1)
			{
				z$phase1Its <- z$phase1Its + 1
			}
		}
		else ### only do parallels at this level for forward simulation
		{
			zz <- clusterCall(z$cl, simstats0c, zsmall, xsmall)
			z$n <- z$n + z$int
			if (phase == 1)
			{
				z$phase1Its <- z$phase1Its + z$int
			}
		}
		## now update the stores
		if (z$int == 1)
		{
			fra <- colSums(zz$fra)
			fra <- fra - z$targets
			fra2 <- zz$fra
			z$sf[z$nit, ] <- fra
			z$sf2[z$nit, , ] <- zz$fra
			z$sims[[z$nit]] <- zz$sims
			z$chain[[z$nit]] <- zz$chain
			fra <- fra + z$targets
		}
		else
		{
			for (i in 1:int)
			{
				fra <- colSums(zz[[i]]$fra)
				fra <- fra - z$targets
				z$sf[z$nit + (i - 1), ] <- fra
				z$sf2[z$nit + (i - 1), , ] <- zz[[i]]$fra
				z$sims[[z$nit + (i - 1)]] <- zz[[i]]$sims
			}
			fra2 <- t(sapply(zz, function(x)x$fra))
			dim(fra2) <- c(int, nrow(zz[[1]]$fra), z$pp)
			fra <- t(sapply(zz, function(x) colSums(x$fra)))
		}
		if (x$maxlike)
		{
			z$sdf[[z$nit]] <- zz$dff
			z$sdf2[[z$nit]] <- zz$dff2
			z$accepts[z$nit, , ] <- zz$accepts
			z$rejects[z$nit, , ] <- zz$rejects
			z$aborts[z$nit, , ] <- zz$aborts
		}
		else if (z$FinDiff.method)
		{
			z <- FiniteDifferences(z, x, fra, fra2)
			for (i in 0:(z$int - 1))
			{
				z$sdf[[z$nit + i]] <- z$sdf0[i + 1, , ]
				if (z$byWave)
				{
					z$sdf2[[z$nit + i]] <- z$sdf02[i + 1, , , ]
				}
			}
		}
		else
		{
			if (int==1)
			{
				if (!is.null(zz[['sc']]))
				{
					z$ssc[z$nit , ,] <- zz$sc
				}
			}
			else
			{
                for (i in 1:int)
                {
                    if (!is.null(zz[[i]][['sc']]))
					{
                        z$ssc[z$nit + (i - 1), , ] <- zz[[i]]$sc
					}
                }
            }
		}
		if ((!x$maxlike) && z$cconditional && z$Phase == 3)
		{
			if (int == 1)
			{
				z$ntim[z$nit, ] <- zz$ntim0
			}
			else
			{
				for (i in 1:int)
				{
					z$ntim[z$nit + (i-1), ] <- zz[[i]]$ntim0
				}
			}
		}
		## post processing also differs
		if (phase == 1)
		{
			CheckBreaks()
			if (UserInterruptFlag() || UserRestartFlag())
			{
				return(z)
			}
			val <- getProgressBar(z$pb)
			if (z$FinDiff.method)
			{
				val <- val + z$pp + 1
			}
			else
			{
				val <- val + 1
			}
			z$pb <- setProgressBar(z$pb, val)
			progress <- val / z$pb$pbmax * 100
			if (z$nit <= 5 || z$nit %% z$writefreq == 0 || z$nit %%5 == 0 ||
				x$maxlike || z$FinDiff.method ||
				(z$int > 1 && z$nit %% z$writefreq < z$int))
			{
				##  Report(c('Phase', z$Phase, 'Iteration ', z$nit, '\n'))
				if (is.batch())
				{
					Report(c('Phase ', z$Phase, ' Iteration ', z$nit,
							 ' Progress: ', round(progress), '%\n'), sep='')
				}
				else
				{
					DisplayTheta(z)
					DisplayDeviations(z, z$sf[z$nit, ])
				}
			}
			if (!z$OK || UserRestartFlag() || UserInterruptFlag())
				return(z)
		}
		else ## phase 3
		{
			if (!is.batch())
			{
				if (nit < 10)
				{
					Report(c("  ", nit, " ",
							 format(z$sf[nit,], width=1, digits=4, nsmall=4),
							 "\n"), cf)
				}
				if (nit >= 10)
				{
					CheckBreaks()
					##    if (nit==10) set up stopkey hint
					##   Early termination of estimation
					if( UserInterruptFlag())
					{
						Report(c("The user asked for an early stop of the ",
								 "algorithm during phase 3.\n"), outf)
						z$Phase3Interrupt <- TRUE
						if (nit < 500)
						{
							if (EarlyEndPhase2Flag())
							{
								Report('This implies that ', outf)
							}
							else
							{
								Report(c("This implies that the estimates",
										 "are as usual,\nbut the "), outf)
							}
							Report(c("diagnostic checks, covariance",
									 "matrices and t-values \nare less ",
									 "reliable, because they ",
									 'are now based on only', nit,
									 'phase-3 iterations.\n\n'), outf)
						}
						z$sf <- z$sf[1:nit, , drop=FALSE]
						z$sf2 <- z$sf2[1:nit, , , drop=FALSE]
						if (!x$maxlike)
						{
							z$ssc <- z$ssc[1:nit, , , drop=FALSE]
						}
						else
						{
							z$sdf <-z$sdf[1:nit]
							z$sdf2 <-z$sdf2[1:nit]
						}
						z$sims <-z$sims[1:nit]
						z$Phase3nits <- nit
						z$n3 <- nit
						break
					}
					if (UserRestartFlag())
						break
				}
			}
		}
		if (!z$OK || UserRestartFlag())
		{
			return(z)
		}
	}
	z
}

