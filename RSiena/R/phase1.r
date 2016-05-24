##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: phase1.r
## *
## * Description: Phase 1 processing is controlled by the functions phase1.1
## * and phase1.2. phase1.1 does the first 10 iterations and checks that
## * finite differencing is going OK. If not, alters epsilon and forces a
## * restart. Phase 1.2 does the rest of the iterations and then calculates
## * the derivative estimate. Both call doPhase1or3iterations which does a block
## * of iterations, to save copying large objects around.
## ****************************************************************************/
##args: x model object (readonly), z control object
##
##@phase1.1 siena07 Do first 10 iters (before check if using finite differences)
phase1.1 <- function(z, x, ...)
{
    ## initialise phase 1
    z$SomeFixed <- FALSE
    z$Phase <-  1
    f <- FRANstore()
    int <- z$int
    z <- AnnouncePhase(z, x)
    z$phase1Its <- 0
    ## fix up iteration numbers if using multiple processors
    if (10 %% int == 0)
    {
        firstNit <- 10
    }
    else
    {
        firstNit <- 10 + int - 10 %% int
    }
    if ((z$n1 - firstNit) %% int == 0)
    {
        endNit <- z$n1
    }
    else
    {
        endNit <- z$n1  + int - (z$n1 - firstNit) %% int
    }
    z$n1 <- endNit
	z <- createSiena07stores(z, z$n1, f)
    z$writefreq <- 10
    z$DerivativeProblem <- FALSE
    z$Deriv <- !z$FinDiff.method ## can do both in phase 3 but not here!
    xsmall <- NULL
    zsmall <- makeZsmall(z)

    nits <- seq(1, firstNit, int)
    if (any(nits >= 6))
    {
        nits6 <- min(nits[nits >= 6 ])
    }
    else
	{
        nits6 <- -1
	}
	### call the subroutine to do the iterations
	z <- doPhase1or3Iterations(1, z, x, zsmall, xsmall, nits, nits6)

	if (UserInterruptFlag() || UserRestartFlag())
	{
		return(z)
	}
    if (z$FinDiff.method)
    {
        npos <- z$npos
        if (any(npos[!z$fixed] < 5))
        {
            j<- (1 : z$pp)[!z$fixed & npos < 5]
            for (i in 1 : length(j))
            {
                Report(c("After ", z$n, " iterations, only ", npos[j[i]],
                         " differences in coordinate ", j[i], ".\n"), sep="",
                       cf)
                Report(c("with epsilon = ", z$epsilon[j[i]], ".\n"), sep="", cf)
            }
            use <- !z$fixed & npos <= 2
            z$epsilon[use] <- ifelse(z$posj[use],
                                    3.0 * z$epsilon[use],
                                    10.0 * z$epsilon[use])
            use <- !z$fixed & npos > 2 & npos < 5
            z$epsilon[use] <- ifelse(z$posj[use],
                                  2.0 * z$epsilon[use],
                                  sqrt(10.0) * z$epsilon[use])
            use <- !z$fixed & npos < 5
            z$epsilon[use] <- pmin(100.0 * z$scale[use], z$epsilon[use])
            z$epsilon[use] <- pmax(0.1 * z$scale[use], z$epsilon[use])
            Report(c("New epsilon =", paste(" ", z$epsilon[use], collapse=""),
					 ".\n"), sep="", cf, fill=80)
            if (z$repeatsForEpsilon <= 4)
            {
                Report("Change value of epsilon and restart Phase 1.\n", cf)
                z$epsilonProblem <- TRUE
                return(z)
            }
           if (any(npos[!z$fixed] <= 1))
            {
                j0s <- which(!z$fixed & npos <= 1)
                for (j0 in j0s)
                {
                    Report (c("Difficulties with parameter", j0, ".\n"),
                            outf)
                    Report(c("After ", z$n, " iterations, only", npos[j0],
                             "differences in this coordinate.\n"), outf)
                    Report("This parameter is fixed and not estimated.\n",
                           outf)
                    Report(c("Fix parameter", j0, ", in phase 1.\n"), cf)
                }
                z$newFixed[j0s] <- TRUE
                z$fixed[j0s] <- TRUE
                z$SomeFixed <- TRUE
            }
            if (all(z$fixed))
            {
                z$AllNowFixed <- TRUE
            }
        }
    }
    z
}


##@phase1.2 siena07 Do rest of phase 1 iterations
phase1.2 <- function(z, x, ...)
{
    ##finish phase 1 iterations and do end-of-phase processing
    zsmall <- makeZsmall(z)
    xsmall <- NULL
    int <- z$int
    if (z$n1 > z$phase1Its)
    {
        nits <- seq((z$phase1Its+1), z$n1, int)
		z <- doPhase1or3Iterations(1, z, x, zsmall, xsmall, nits)
    }
	if (UserInterrupt() || UserRestartFlag())
	{
		return(z)
	}
	# for dolby option: regress deviations fra on scores ssc
	z$regrCoef <- rep(0, z$pp)
	z$regrCor <- rep(0, z$pp)
	if (!is.null(z$ssc))
	{
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
		Report('Regression coefficients:\n', cf)
		PrtOutMat(format(as.matrix(t(z$regrCoef)), digits = 2, nsmall = 2), cf)
	}

    z$timePhase1 <- (proc.time()['elapsed'] - z$ctime) / (z$nit - 1)
    if (x$checktime  && !is.na(z$timePhase1))
    {
        Report(c('Time per iteration in phase 1  =',
                 format(z$timePhase1, digits = 4, nsmall = 4),'\n'), lf)
    }
    z$mnfra <- colMeans(z$sf)
    Report('Average deviations NR generated statistics and targets\n', cf)
    Report('after phase 1:\n', cf)
    PrtOutMat(format(as.matrix(z$mnfra), width = 15, nsmall = 6), cf)
    z <- CalculateDerivative(z, x)
	## browser()
    if (!z$OK || z$DerivativeProblem) ##longer phase 1 or use finite differences
    {
        return(z)
    }
    z$SomeFixed <-z$SomeFixed | z$cdSomeFixed
    if (z$SomeFixed)
    {
        Report('(Values for fixed parameters are meaningless.)\n', cf)
        z$dfra[outer(z$fixed, z$fixed, '|')] <- 0
        diag(z$dfra)[z$fixed] <- 1.0
        z$mnfra[z$fixed] <- 0.0
        z$sf[ , z$fixed] <- 0
    }
	# Manage derivative matrix
    if (any(diag(z$dfra) < 1e-8))
    {
        sml <- (1 : z$pp)[diag(z$dfra) < 1e-8]
        diag(z$dfra)[sml] <- 1e-3
        Report(c('Diagonal elements(s)', sml,
                 'of derivative matrix amended in phase 1.'), cf, fill=80)
    }
	# Invert derivative matrix and define step fchange that will not be made.
    if (inherits(try(dinv <- solve(z$dfra)), "try-error"))
    {
        Report('Error message for inversion of dfra: \n', cf)
        diag(z$dfra) <- diag(z$dfra) + 1
        Report('Intervention 1.4.1: ridge added after phase 1.\n', cf)
        if (inherits(try(dinv <- solve(z$dfra)), "try-error"))
        {
            Report(c('Error. After phase 1, derivative matrix non-invertible',
                     'even with a ridge.\n'), cf)
            fchange <- 0
			stop("Cannot proceed: derivative matrix not invertible")
        }
        else
        {
            fchange <- as.vector(dinv %*% z$mnfra)
            z$dinv <- dinv
        }
    }
    else
    {
        fchange <- as.vector(dinv %*% z$mnfra)
        z$dinv <- dinv
    }
	if (!x$diagg) 
	{
		# Partial diagonalization of derivative matrix
		# for use if 0 < x$diagonalize < 1.
		temp <- (1-x$diagonalize)*z$dfra + x$diagonalize*diag(diag(z$dfra))
		temp[z$fixed, ] <- 0.0
		temp[, z$fixed] <- 0.0
		diag(temp)[z$fixed] <- 1.0
		# Invert this matrix
		z$dinvv <- solve(temp)
	}
    Report('dfra :\n', cf)
    ##  browser()
    PrtOutMat(z$dfra, cf)
    Report('inverse of dfra :\n', cf)
    PrtOutMat(dinv, cf)
    Report('dinvv :\n', cf)
    PrtOutMat(z$dinvv, cf)
    Report('Full Quasi-Newton-Raphson step after phase 1:\n', cf)
    Report(c(paste(1:z$pp, '. ', format(-fchange, width = 12, digits = 6,
                                        nsmall = 6), sep = '', collapse = '\n'),
             '\n'), cf)
    Report(c('This step is multiplied by the factor ',
             format(0.5 * x$firstg, digits = 5, nsmall = 5, width = 8), '.\n'),
           cf, sep='')
    fchange <- 0.5 * x$firstg * fchange
    fchange[z$fixed] <- 0.0
    ##check if jump is too large
    maxrat<- max(abs(fchange / z$scale))
    if (maxrat > 10.0)
    {
        fchange <- 10 * fchange / maxrat;
        Report(c('Intervention 1.4.2: jump after phase 1 decreased by factor',
                 maxrat, '.\n'), cf)
    }
    ##check positivity
    if (z$anyposj)
    {
        neg <- z$posj & fchange >= z$theta
        if (any(neg))
        {
            fchange[neg] <- 0.5 * z$theta[neg]
            Report(c('Intervention 1.4.3: positivity restriction after phase 1',
                     'coordinate(s)', (1 : z$pp)[neg], '.'), cf, fill=80)
            Report('\n', cf)
        }
    }
    ## make step
    if (x$nsub >= 1)
    {
        z$theta[!z$fixed] <- z$theta[!z$fixed] - fchange[!z$fixed]
    }
    Report(c('Phase 1 achieved after ', z$phase1Its, ' iterations.\n'), cf)
    WriteOutTheta(z)
    z$nitPhase1 <- z$phase1Its
    z$phase1devs <- z$sf
    z$phase1dfra <- z$frda
    z$phase1sdf <- z$sdf
    z$phase1sdf2 <- z$sdf2
    z$phase1scores <- z$ssc
    z$phase1accepts <- z$accepts
    z$phase1rejects <- z$rejects
	z$phase1aborts <- z$aborts
    ##browser()
    z
}

##@CalculateDerivative siena07 Calculates derivative in Phase 1
CalculateDerivative <- function(z, x)
{
    if (z$FinDiff.method || x$maxlike)
    {
  		dfra <- t(as.matrix(Reduce("+", z$sdf) / length(z$sdf)))
	}
    else
    {
        ##note that warnings is never set as this piece of code is not executed
        ##if force_fin_diff_phase_1 is true so warning is zero on entry
        ##browser()

        dfra <- derivativeFromScoresAndDeviations(z$ssc, z$sf2)

        z$jacobianwarn1 <- rep(FALSE, z$pp)
        if (any(diag(dfra) <= 0))
         {
            for (i in 1 : z$pp)
            {
                if (dfra[i, i] < 0)
                {
                    ##browser()
                    z$jacobianwarn1[i] <- TRUE
                    Report(c('Warning: diagonal value', i,
                             'is non-positive.\n\n'), cf)
                }
            }
            if (z$n1 < 200)
            {
                z$n1 <- z$n1 * 2
                Report(c('New phase 1 of increased length because of',
                         'unreliable derivative matrix\n'), cf)
                Report(c('New phase 1 of increased length because of',
                         'unreliable derivative matrix\n'), lf)
                z$DerivativeProblem <- TRUE
                z$LongerPhase1 <-  TRUE
            }
            else
            {
                z$FinDiff.method <- TRUE
                Report(c('New phase 1 with finite differences because of',
                         'unreliable derivative matrix\n'), cf)
                Report(c('New phase 1 with finite differences because of',
                         'unreliable derivative matrix\n'), lf)
               z$DerivativeProblem <- TRUE
            }
            return(z)
        }
    }
    dfra[outer(z$fixed,z$fixed,'|')] <- 0
    diag(dfra)[z$fixed] <- 1.0
    Report('Diagonal values of derivative matrix :\n', cf)
    Report(format(diag(dfra), digits = 4, nsmall = 4, width = 8), cf, fill=80)
    someFixed <- FALSE
    if (any(diag(dfra) <= 0 & !z$fixed))
    {
        neg<- which(diag(dfra) <= 0 & !z$fixed)
        z$fixed[neg] <- TRUE
        z$newFixed[neg] <- TRUE
        someFixed <- TRUE
        Report(c('*** nonpositive diagonal value(s):', neg,
                  ' is/are fixed.\n'), cf)
            Report(c('Estimation problem with parameter(s)', neg,
                     'this/these parameter(s) is/are fixed.\n'), outf)
    }
    z$dfra <- dfra
    z$cdSomeFixed <- someFixed
    z
}

##@FiniteDifferences siena07 Does the extra iterations for finite differences
FiniteDifferences <- function(z, x, fra, fra2)
{
    int <- z$int
    fras <- array(0, dim = c(int, z$pp, z$pp))
	if (z$byWave)
	{
		fras2 <- array(0, dim = c(int, z$f$observations - 1, z$pp, z$pp))
	}
	else
	{
		fras2 <- NULL
	}
    xsmall <- NULL
    for (i in 1 : z$pp)
    {
        zdummy <- makeZsmall(z)
        if (z$Phase == 3 || !z$fixed[i])
        {
            zdummy$theta[i] <- z$theta[i] + z$epsilon[i]
        }
        if (int == 1)
        {
            zz <- x$FRAN(zdummy, xsmall, fromFiniteDiff=TRUE)
            if (!zz$OK)
            {
                z$OK <- zz$OK
                return(z)
            }
        }
        else
        {
            zz <- clusterCall(z$cl, simstats0c, zdummy, xsmall,
                              fromFiniteDiff=TRUE)
        }
        if (int == 1)
        {
            fras[1, i, ] <- colSums(zz$fra) - fra
			if (z$byWave)
			{
				fras2[1, , i, ] <- zz$fra - fra2
			}
        }
        else
        {
            for (j in 1 : int)
            {
                fras[j, i, ] <- colSums(zz[[j]]$fra) - fra[j, ]
				if (z$byWave)
				{
					fras2[j, , i, ] <- zz[[j]]$fra - fra2[j, , ]
				}
           }
        }
    }
    if (z$Phase == 1 && z$nit <= 10)
    {
        for (ii in 1: min(10 - z$nit + 1, int))
        {
            z$npos <- z$npos +
                ifelse(abs(diag(matrix(fras[ii, , ], nrow=z$pp))) > 1e-6, 1, 0)
        }
    }
    z$sdf0 <- fras / rep(rep(z$epsilon, each=int), z$pp)
	if (z$byWave)
	{
		z$sdf02 <- fras2 / rep(rep(z$epsilon, each=int * dim(fras2)[2]), z$pp)
	}
    z
}
##@derivativeFromScoresAndDeviations siena07 create dfra from scores and deviations
derivativeFromScoresAndDeviations <- function(scores, deviations)
{
    nIterations <- dim(scores)[1]
    nWaves <- dim(scores)[2]
    nParameters <- dim(scores)[3]
    ## replaced nested applies by memory saving for loops after problems
    dfra <- matrix(0, nrow=nParameters, ncol=nParameters)
    for (i in 1:nIterations)
    {
        for (j in 1:nWaves)
        {
            dfra <- dfra + outer(scores[i, j, ], deviations[i, j, ])
        }
    }
    dfra <- t(dfra)
    dfra<- dfra / nIterations
    tmp <- matrix(sapply(1 : nWaves, function(i)
                         outer(colMeans(deviations)[i,],
                               colMeans(scores)[i,])), ncol=nWaves)

    dfra - matrix(rowSums(tmp), nrow=nParameters)
}

##@makeZsmall siena07 create a minimal version of z to pass between processors.
makeZsmall <- function(z)
{
    zsmall <- NULL
    zsmall$theta <- z$theta
    zsmall$Deriv <- z$Deriv
    zsmall$Phase <- z$Phase
    zsmall$nit <- z$nit
    zsmall$FinDiff.method <- z$FinDiff.method
    zsmall$int2 <- z$int2
    zsmall$cl <- z$cl
    zsmall$maxlike <- z$maxlike
    zsmall$cconditional <- z$cconditional
    zsmall$condvar <- z$condvar
    zsmall$pp <- z$pp
    zsmall$nrunMH <- z$nrunMH
    zsmall$returnDeps <- z$returnDeps
    zsmall$returnChains <- z$returnChains
    zsmall$returnDataFrame <- z$returnDataFrame
    zsmall$addChainToStore <- z$addChainToStore
	zsmall$callGrid <- z$callGrid
	zsmall$thetaMat <- z$thetaMat
	zsmall$byWave <- z$byWave
    zsmall
}

##@createSiena07Stores siena07 set up the storage areas used in phase 1 and 3
createSiena07stores <- function(z, nIterations, f)
{
    z$sf <- matrix(0, nrow = nIterations , ncol = z$pp)
    z$sf2 <- array(0, dim=c(nIterations, f$observations - 1, z$pp))
	if (!z$maxlike && !z$FinDiff.method)
	{
		z$ssc <- array(0, dim=c(nIterations, f$observations - 1, z$pp))
	}
	else
	{
		z$sdf <- vector("list", nIterations)
		z$sdf2 <- vector("list", nIterations)
	}
	## misdat steps are separated out giving 9 types
    z$accepts <- array(0, dim=c(nIterations, z$nDependentVariables, 9))
    z$rejects <- array(0, dim=c(nIterations, z$nDependentVariables, 9))
    z$aborts <- array(0, dim=c(nIterations, z$nDependentVariables, 9))
	z$npos <- rep(0, z$pp)
    if (!is.null(z$cconditional) && z$cconditional)
    {
        z$ntim <- matrix(NA, nrow=nIterations, ncol=f$observations - 1)
    }
    z$sims <- vector("list", nIterations)
	z
}
