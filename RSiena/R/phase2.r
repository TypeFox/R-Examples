##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: phase2.r
## *
## * Description: This module contains the functions phase2.1, proc2subphase
## * and doIterations which together perform a robbins-monro stochastic
## * approximation algorithm.
## * phase2.1 and proc2subphase are called from robmon in robmon.r.
## ****************************************************************************/
## args: z: internal control object
##       x: model object (readonly as not returned)

##@storeinFRANstore siena07 Used to avoid Namespace problems with multiple processes
storeinFRANstore <- function(...)
{
    FRANstore(...)
}
##@phase2.1 siena07 Start phase 2
phase2.1<- function(z, x, ...)
{
    #initialise phase2
    if (x$maxlike)
    {
        z$phase2fras <- array(0, dim=c(4, z$pp, 1000))
        z$rejectprops <- array(0, dim=c(4, 4, 1000))
    }
    z$Phase <- 2
    z$writefreq <- 1
    if (!is.batch())
    {
        tkconfigure(z$tkvars$earlyEndPhase2,state='normal')
        tkconfigure(z$tkvars$subphase,state='normal')
        tkconfigure(z$tkvars$subphaselabel,state='normal')
    }
    z$Deriv <- FALSE
    z$sd <- sqrt(apply(z$sf, 2, function(x) sum(x^2) / nrow(z$sf) - mean(x)^2))
    z$sd[z$fixed] <- 0
    Report(paste('\nPhase 2 has', x$nsub, 'subphases.\n'), cf)
    z$gain <- x$firstg
    if (x$nsub <= 0)
    {
        Report('With 0 subphases, there is no phase 2.\n', cf)
    }
    else
    {
        if (z$maxrepeatsubphase >= 2)
        {
            Report(c('Each subphase can be repeated up to', z$maxrepeatsubphase,
                   'times\n'), cf)
        }
        z <- proc2subphase(z, x, 1, ...)
    }
    z
}
##@proc2subphase siena07 Do one subphase of phase 2
proc2subphase <- function(z, x, subphase, useAverage=TRUE, ...)
{
    ## init subphase of phase 2
    z <- AnnouncePhase(z, x, subphase)
	Report(paste("\nStart phase ", z$Phase, ".", subphase, "\n", sep=""), cf)
	if (subphase <= 0)
	{
		z$n2min <- 5
		z$n2max <- 100
	}
	else
	{
		z$n2min <- z$n2minimum[subphase]
		z$n2max <- z$n2maximum[subphase]
	}
    z$repeatsubphase <- 0
    repeat
    {
        z$repeatsubphase <- z$repeatsubphase + 1
        z$truncated <- rep(FALSE, z$n2max)
        z$positivized <- matrix(FALSE, nrow = z$n2max, ncol = z$pp)
        z$ctime <- proc.time()[3]
        z$time1 <- proc.time()[3]
        z$thav <- z$theta
        z$thavn <- 1
        ## cat(z$thav, z$theta, '\n')
        z$prod0 <- rep(0, z$pp)
        z$prod1 <- rep(0, z$pp)
        ## ###############################################
        ## do the iterations for this repeat of this subphase
        ## ##############################################
        ##z <- doIterationsCopy(z, x, subphase, ...) removed as out of sync
        z <- doIterations(z, x, subphase, ...)
        ##   if (z$nit == 50) browser()
        if (!z$OK || UserInterruptFlag() || UserRestartFlag() ||
            EarlyEndPhase2Flag())
        {
            break
        }
        ##
        ## end processing for this repeat of this subphase
        ##
        ## report truncations and positivizations
        if (any(z$truncated))
        {
            msg<- paste('Intervention 2.', subphase,
                        '.1: changes truncated, iterations: ',sep='')
            Report(msg, cf)
            Report(which(z$truncated), cf, fill=80)
        }
        if (any(z$positivized))
        {
            msg <- paste('Intervention 2.',subphase,
                         '.2: positivity restriction:\n ',sep='')
            Report(msg,cf)
            subs <- which(rowSums(z$positivized) > 0)
            msg<- sapply(subs, function(i, y)
                         paste('Observation:', i, 'Coordinate(s):',
                               paste((1:z$pp)[y[i,]], collapse = ' ')),
                         y = z$positivized)
            Report(msg, cf, fill=80)
        }
        if ((subphase >= 1) && (z$maxacor >= sqrt(2.0 / (z$nit + 1))))
		{
			Report('Warning: an autocorrelation is positive at the end',cf)
			Report(' of this subphase.\n',cf)
			Report('Autocorrelations:\n',cf)
			prtmat <- z$prod1 / z$prod0
			PrtOutMat(as.matrix(prtmat), cf)
        }
        if ((z$nit >= z$n2max)
			|| (subphase <= 0)
			|| (z$maxacor < 1e-10)
			|| (z$repeatsubphase >= z$maxrepeatsubphase))
        {
            break
        }
    }
   ##finalize the subphase
    if (!z$OK || UserInterruptFlag() || UserRestartFlag())
    {
        return(z)
    }
    if (EarlyEndPhase2Flag())
    {
        Report('The user asked for early end of phase 2.\n', outf)
    }
    if (useAverage)
    {
		z$theta <- z$thav / z$thavn   # z$thavn = (z$nit + 1)
    }
    else
	{
		cat('\n')
		cat('Regression\n')
		##  use regression
		cat(z$thav / z$thavn, '; ') #(z$nit + 1)
		mylm <- lm(z$sf[1:z$nit, ] ~ z$thetaStore[1:z$nit, ])
		coefs <- coef(mylm)[-1, ]
		newvals <- solve(t(coefs), - mylm$coef[1, ])
		z$theta <- newvals
		cat(z$theta, '\n')
	}
    DisplayThetaAutocor(z)
    ##    cat('it',z$nit,'\n')
    ##recalculate autocor using -1 instead of -2 as error
    ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -1)
    maxacor <- max(-99, ac[!z$fixed]) ##note -1 > -99
    Report(paste('Phase ', z$Phase,'.', subphase, ' ended after ', z$nit,
                 ' iterations.\n', sep = ''), cf)
    if (x$checktime)
    {
        time1 <- proc.time()[[3]] - z$ctime
        subphaseTime <- time1 / z$nit
        Report(paste('Time per iteration in phase ', z$Phase, '.', subphase,
                     ' = ', format(subphaseTime, nsmall=4, digits=4),
                     '\n', sep=''), lf)
    }
    if ((maxacor >= sqrt(2 / (z$nit + 1))) && (subphase >= 1))
    {
        Report('Warning. Autocorrelation criterion not satisfied.\n', cf)
    }
    WriteOutTheta(z)
    if (EarlyEndPhase2Flag())
    {
        return(z)
    }
    if (subphase == 2 && z$restart) ## this means we restarted in phase 1
        ##because of epsilon and need to restart the whole thing again now
    {
        Report('Restart after subphase 2.2 from current parameter values\n', cf)
        Report('because the initial values used ', cf)
        Report('led to questionable epsilon values\n', cf)
        z$fixed[z$newfixed] <- FALSE
        z$restarted <- TRUE
    }
	## For the next subphase:
	if (x$maxlike)
	{
		z$gain <- z$gain * 0.25
	}
	else
	{
		z$gain <- z$gain * 0.5
	}
    z
} ##end of this subphase

##@doIterations siena07 Do all iterations for 1 repeat of 1 subphase of phase 2
doIterations<- function(z, x, subphase,...)
{
    z$nit <- 0
    ac <- 0
	z$maxacor <- 1
    xsmall <- NULL
    zsmall <- makeZsmall(z)
    z$returnDeps <- FALSE
    repeat
    {
        z$n <- z$n+1
        z$nit <- z$nit + 1
        if (subphase == 1 && z$nit == 2)
            z$time1 <- proc.time()[[3]]
        if (subphase == 1 && z$nit == 11)
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
            z$writefreq <- roundfreq(z$writefreq)
        }
        if ((z$nit <= 10) || (z$nit %% z$writefreq ==0))
        {
            DisplayIteration(z)
           if (is.batch())
            {
                val <- getProgressBar(z$pb)
                increment <- ifelse(z$nit <= 10, 1, z$writefreq)
                Report(paste('Phase ', z$Phase, ' Subphase ', subphase,
                             ' Iteration ', z$nit,' Progress: ',
                             round((increment + val) / z$pb$pbmax * 100),
                             '%\n', sep = ''))
                z$pb <- setProgressBar(z$pb, val + increment)
            }
            else
            {
               if (z$nit>1)
                {
                    DisplayDeviations(z, fra)
                }
                if  (z$nit %% z$writefreq == 0)
                {
                    val <- getProgressBar(z$pb)
                    z$pb <-setProgressBar(z$pb, val + z$writefreq)
                }
            }
        }
        zsmall$nit <- z$nit
		if (x$dolby) {zsmall$Deriv <- TRUE} ## include scores in FRAN
        if (z$int == 1) ## then no parallel runs at this level
        {
            zz <- x$FRAN(zsmall, xsmall)
            fra <- colSums(zz$fra) - z$targets
            if (!zz$OK)
            {
                z$OK <- zz$OK
                break
            }
			if (x$dolby)
			## subtract regression on scores;
			## permitted because fra is used in phase2 only for updates
			{
			ssc <- zz$sc # scores;  periods by variables
			ssc <- colSums(ssc) # add over periods
			fra <- fra - (z$regrCoef * ssc)
			}
        }
        else
        {
            zz <- clusterCall(z$cl, simstats0c, zsmall, xsmall)
			if (x$dolby)
			## subtract regression on scores;
			{
				fra <- sapply(zz, function(x){
					ssc <- x$sc # scores;  periods by variables
					ssc <- colSums(ssc) # add over periods
					colSums(x$fra)- z$targets - (z$regrCoef * ssc)
					})
			}
			else
			{
				fra <- sapply(zz, function(x) colSums(x$fra)- z$targets)
			}
            dim(fra) <- c(z$pp, z$int)
            fra <- rowMeans(fra)
            zz$OK <- sapply(zz, function(x) x$OK)
            if (!all(zz$OK))
            {
                z$OK <- FALSE
                break
            }

        }
		## setup check for end of phase. either store or calculate
        if (z$nit %% 2 == 1)
        {
            prev.fra <- fra
        }
        else
        {
            z$prod0 <- z$prod0 + fra * fra
            z$prod1 <- z$prod1 + fra * prev.fra
            ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -2)
            z$maxacor <- max(-99, ac[!z$fixed]) ##note -2 > -99
            z$minacor <- min(1, ac[(!z$fixed) & ac > -1.0])
            z$ac <- ac
            if  (z$nit %% z$writefreq == 0)
            {
                DisplayThetaAutocor(z)
            }
        }
        ## limit change.  Reporting is delayed to end of phase.
# The old version had a numerical ifelse condition:
#          maxrat<- max(ifelse(z$sd, abs(fra)/ z$sd, 1.0))
# But ifelse with a numerical first argument only tests the value 0,
# which here has probability 0... So this really has no effect.
# I (TS) do not understand this.
# Therefore I changed it.
        if (x$diagg)
		{
            changestep <- fra / diag(z$dfra)
		}
		else
		{
			changestep <- as.vector(fra %*% z$dinvv)
		}
		changestep[z$fixed] <- 0.0
		maxratt <- max(abs(changestep))
		if (maxratt > 5)
		{
			changestep <- 5*changestep/maxratt
           	z$truncated[z$nit] <- TRUE
		}
        fchange <- as.vector(z$gain * changestep)
        ## check positivity restriction
        z$positivized[z$nit, ] <- z$posj & (fchange > z$theta)
        fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)
        zsmall$theta <- zsmall$theta - fchange
        z$theta <- zsmall$theta
#		z$thetaStore[z$nit,] <- z$theta
        z$thav <- z$thav + zsmall$theta
        z$thavn <- z$thavn + 1
        if (x$maxlike && !is.null(x$moreUpdates) && x$moreUpdates > 0)
        {
            z <- doMoreUpdates(z, x, x$moreUpdates * subphase)
            zsmall$theta <- z$theta
        }
        ##check for user interrupt
        CheckBreaks()
        if (UserInterruptFlag() || UserRestartFlag() || EarlyEndPhase2Flag())
        {
            break
        }
        ## do we stop?
        if ((z$nit >= z$n2min && z$maxacor < 1e-10)||
			   (z$nit >= z$n2max)||
			   ((z$nit >= 50) && (z$minacor < -0.8) &&
                (z$repeatsubphase < z$maxrepeatsubphase)))
        {
            break
        }
    }
    z
}

