### function defining the EM algorithm for the CFR estimates

##' A function to estimate the relative case fatality ratio when reporting rates
##' are time-varying and deaths are lagged because of survival time.
##' 
##' This function implements an EM algorithm to estimate the relative case
##' fatality ratio between two groups when reporting rates are time-varying and
##' deaths are lagged because of survival time.
##' @usage EMforCFR(assumed.nu, alpha.start.values, full.data, max.iter = 50, 
##'   verb = FALSE, tol = 1e-10, SEM.var = TRUE)
##'   
##' @param assumed.nu  a vector of probabilities corresponding to the survival
##'   distribution, i.e. nu[i]=Pr(surviving i days | fatal case)
##' @param alpha.start.values a vector starting values for the reporting rate
##'   parameter of the GLM model. This must have length which corresponds to one
##'   less than the number of unique integer values of full.dat[,"new.times"].
##' @param full.data  A matrix of observed data. See description below.
##' @param max.iter The maximum number of iterations for the EM algorithm and
##'   the accompanying SEM algorithm (if used).
##' @param verb An indicator for whether the function should print results as it
##'   runs.
##' @param tol A tolerance to use to test for convergence of the EM algorithm.
##' @param SEM.var If TRUE, the SEM algorithm will be run in addition to the EM
##'   algorithm to calculate the variance of the parameter estimates.
##'   
##' @details The data matrix full.data must have the following columns: 
##'   \describe{ \item{grp}{a 1 or a 2 indicating which of the two groups, j, 
##'   the observation is for.} \item{new.times}{an integer value representing
##'   the time, t, of observation.} \item{R}{the count of recovered cases with
##'   onset at time t in group j.} \item{D}{the count of deaths which occured at
##'   time t in groupo j (note that these deaths did not have disease onset at
##'   time t but rather died at time t).} \item{N}{the total cases at t, j, or
##'   the sum of R and D columns.} }
##'   
##' @return A list with the following elements \describe{ \item{naive.rel.cfr
##'   }{the naive estimate of the relative case fatality ratio} 
##'   \item{glm.rel.cfr }{the reporting-rate-adjusted estimate of the relative
##'   case fatality ratio} \item{EM.rel.cfr }{the lag-adjusted estimate of the
##'   relative case fatality ratio} \item{EM.re.cfr.var }{the variance for the
##'   log-scale lag-adjusted estimator taken from the final M-step} 
##'   \item{EM.rel.cfr.var.SEM }{ the Supplemented EM algorithm variance for the
##'   log-scale lag-adjusted estimator} \item{EM.rel.cfr.chain }{a vector of the
##'   EM algorithm iterates of the lag-adjusted relative CFR estimates} 
##'   \item{EMiter}{the number of iterations needed for the EM algorithm to
##'   converge} \item{EMconv}{indicator for convergence of the EM algorithm.  0
##'   indicates all paramters converged within max.iter iterations.  1 indicates
##'   that the estimate of the relative case fatality ratio converged but other
##'   did not.  2 indicates that the relative case fatality ratio did not
##'   converge.} \item{SEMconv}{indicator for convergence of SEM algorithm. 
##'   Same scheme as EMconv.} \item{ests}{ the coefficient estimates for the
##'   model } \item{ests.chain}{ a matrix with all of the coefficient estimates,
##'   at each EM iteration} \item{DM}{the DM matrix from the SEM algorithm} 
##'   \item{DMiter}{a vector showing how many iterations it took for the
##'   variance component to converge in the SEM algorithm} }
##'  @examples        
##'     ## This is code from the CFR vignette provided in the documentation.
##'        
##' data(simulated.outbreak.deaths)
##' min.cases <- 10 
##' N.1 <- simulated.outbreak.deaths[1:60, "N"] 
##' N.2 <- simulated.outbreak.deaths[61:120, "N"] 
##' first.t <- min(which(N.1 > min.cases & N.2 > min.cases)) 
##' last.t <- max(which(N.1 > min.cases & N.2 > min.cases)) 
##' idx.for.Estep <- first.t:last.t 
##' new.times <- 1:length(idx.for.Estep) 
##' simulated.outbreak.deaths <- cbind(simulated.outbreak.deaths, new.times = NA) 
##' simulated.outbreak.deaths[c(idx.for.Estep, idx.for.Estep + 60), "new.times"] <- rep(new.times, + 2)
##' assumed.nu = c(0, 0.3, 0.4, 0.3)
##' alpha.start <- rep(0, 22)
##'        
##' ## caution! this next line may take several minutes (5-10, depanding on 
##' ##    the speed of your machine) to run.
##' \dontrun{cfr.ests <- EMforCFR(assumed.nu = assumed.nu, 
##'                               alpha.start.values = alpha.start, 
##'                               full.data = simulated.outbreak.deaths, 
##'                               verb = FALSE, 
##'                               SEM.var = TRUE, 
##'                               max.iter = 500, 
##'                               tol = 1e-05)}
##' @keywords coarse data
##' @keywords incomplete data
##' @keywords case fatality ratio
##' @keywords infectious disease
##' @export
EMforCFR <- function(assumed.nu, alpha.start.values, full.data,
		     max.iter=50, verb=FALSE, tol=1e-10, SEM.var=TRUE){
	## full.data is the data output from the observe.epidemic() function
	##    with a column specifying the rows to be used for analysis

	#######################################
	## calculate naive and GLM estimates ##
	#######################################

	D.1 <- sum(full.data[full.data[,"grp"]==1&!is.na(full.data[,"new.times"]),"D"])
	D.2 <- sum(full.data[full.data[,"grp"]==2&!is.na(full.data[,"new.times"]),"D"])
	N.1 <- sum(full.data[full.data[,"grp"]==1&!is.na(full.data[,"new.times"]),"N"])
	N.2 <- sum(full.data[full.data[,"grp"]==2&!is.na(full.data[,"new.times"]),"N"])
	naive.rel.cfr <- sum(D.2)/sum(N.2)/(sum(D.1)/sum(N.1))

	dat <- data.frame(full.data)
	unadj.glm.fit <- glm(dat[,"D"] ~ factor(dat[,"new.times"]) + factor(dat[,"grp"]),
			     offset=log(dat[,"N"]), family=poisson)
	glm.rel.cfr <- exp(coef(unadj.glm.fit))[length(coef(unadj.glm.fit))]

	#########################
	## set up EM algorithm ##
	#########################
	proposed.rel.cfr.chain <- rep(0, max.iter)
	n.params <- max(dat[,"new.times"], na.rm=TRUE)+1
	phi.chain <- matrix(nrow=n.params, ncol=max.iter)
	nlag <- length(assumed.nu)

	## convergence codes:
	## 0 = all parameters converge
	## 1 = rCFR converges, but some other parameters do not
	## 2 = rCFR does not converge
	EMconv <- SEMconv <- 0

	## looping variables
	eps <- 1
	j <- 0
	alpha <- alpha.start.values
	while(eps>tol){
		j <- j+1
		############
		## E STEP ##
		############
		dat <- run.Estep(alpha,full.data=full.data,
				 nlag=nlag, assumed.nu=assumed.nu)

		############
		## M STEP ##
		############
		maxLik <- run.Mstep(dat)
		phi.chain[,j] <- maxLik$phi
		alpha <- phi.chain[2:(n.params-1),j]
		proposed.rel.cfr.chain[j] <- exp(phi.chain[n.params,j])

		## check convergence
		if(j>1)	eps <- max((phi.chain[,j]-phi.chain[,j-1])^2)
		if(j==max.iter) {
			## message("*** WARNING: EM agorithm didn't converge ***")
			if((phi.chain[n.params,j]-phi.chain[n.params,j-1])^2 < eps){
				EMconv <- 1
			} else {
				EMconv <- 2
			}
			break
		}
	}

	var.joint <- maxLik$Var

	phi <- phi.chain[,j]
	proposed.rel.cfr <- proposed.rel.cfr.chain[j]
	EMiter <- j
	if(verb) {
		print(paste("naive estimator =", round(naive.rel.cfr, 3)))
		print(paste("GLM estimate unadjusted for lags =", round(glm.rel.cfr, 3)))
		print(paste(EMiter, "iterations for EM convergence"))
	}

	if(SEM.var){
		DM.out <- SEM.variance(full.data=full.data, dat, phi, max.iter, tol, nlag=nlag, alpha.start.values, assumed.nu=assumed.nu)
		DM <- DM.out$DM
		DMiter <- DM.out$DMiter
		if(length(DM.out$loop.idx)>0) {
			loop.idx <- DM.out$loop.idx
			var.joint <- var.joint[-loop.idx,-loop.idx]
		}
		## set SEM convergence codes
		if(any(DMiter==0)) {
			if(DMiter[length(DMiter)]==0) {
				SEMconv <- 2
			} else { SEMconv <- 1 }

		}

		## calculate CFR variance if it converged
		if(SEMconv<2){
			variance.SEM <- var.joint +
				var.joint %*% DM %*% solve(diag(1, ncol(DM))-DM)
			proposed.rel.cfr.var.SEM <- variance.SEM[nrow(DM), nrow(DM)]
			if(proposed.rel.cfr.var.SEM<0) {
				proposed.rel.cfr.var.SEM <- NA
				SEMconv <- 2
			}
		} else { proposed.rel.cfr.var.SEM <- NA }

		if(verb){
			print("Estimate with SEM")
			print(paste("proposed CFR =", round(proposed.rel.cfr, 3),
				    "; 95% CI (",
				    round(exp(log(proposed.rel.cfr)
					      -2*sqrt(proposed.rel.cfr.var.SEM)),3),
				    ",",
				    round(exp(log(proposed.rel.cfr)
					      +2*sqrt(proposed.rel.cfr.var.SEM)),3),
				    ")"))
			if(any(DM.out$DMiter==0)) {
				print("non-convergent DM indices:")
				print(which(DM.out$DMiter==0))
			}
		}
	} else {
		SEMconv <- NA
		DM <- NULL
		DMiter <- NULL
		proposed.rel.cfr.var.SEM <- NULL
	}

	## ##############
	## OUTPUT LIST ##
	## ##############
	out <- list(naive.rel.cfr=naive.rel.cfr,
		    glm.rel.cfr=glm.rel.cfr,
		    EM.rel.cfr=proposed.rel.cfr,
		    EM.rel.cfr.var=var.joint[nrow(var.joint), nrow(var.joint)],
		    EM.rel.cfr.var.SEM = proposed.rel.cfr.var.SEM,
		    EM.rel.cfr.chain=proposed.rel.cfr.chain,
		    EMiter=EMiter,
		    EMconv=EMconv,
		    SEMconv=SEMconv,
		    ests=phi,
		    ests.chain.EM=phi.chain,
		    DM=DM,
		    DMiter=DMiter)

	if(verb) {
		print("Estimate with just GLM variance")
		print(paste("proposed CFR =", round(proposed.rel.cfr, 3),
			    "; 95% CI (",
			    round(exp(log(proposed.rel.cfr)
				      -2*sqrt(out$EM.rel.cfr.var)),3),
			    ",",
			    round(exp(log(proposed.rel.cfr)
				      +2*sqrt(out$EM.rel.cfr.var)),3),
			    ")"))
	}

	return(out)

}


## This function is meant to be run only through the function EMforCFR() 
##   and is used to calculate the variance via the Supplemented EM 
##   algorithm (see Meng and Rubin, 1991)

## params:
##   \item{full.data}{ A matrix of observed data. See description in EMforCFR helpfile.}
##   \item{dat}{A data frame.}
##   \item{phi}{ A vector of fitted parameters from the final EM iteration.}
##   \item{max.iter}{ The maximum number of iterations for SEM algorithm. }
##   \item{tol}{ A tolerance to use to test for convergence of the EM algorithm. }
##   \item{nlag}{ The number of time units for lagged data.  Corresponds to length(assumed.nu).}  
##   \item{alpha.start.values}{ a vector starting values for the reporting rate parameter of the GLM model. This must have length which corresponds to one less than the number of unique integer values of full.dat[,"new.times"].}
##   \item{assumed.nu}{ a vector of probabilities corresponding to the survival distribution, i.e. nu[i]=Pr(surviving i days | fatal case) }

## returned list:
##    \item{DM}{The estimate of the variance-covariance matrix for the model parameters.  Only converged rows are returned.}
##    \item{DMiter}{A vector whose ith entry is the number of iterations needed for convergence of the ith row of the DM matrix.}
##    \item{loop.idx}{If not NULL, the values correspond to the original indices of DM which have been omitted because of lack of convergence.}
SEM.variance <- function(full.data, dat, phi, max.iter, tol, nlag,
			 alpha.start.values, assumed.nu){
	## algorithm parameters
	eps <- 1
	iter <- 0
	phi.hat <- phi ## using phi from notation in Bayesian Data Analysis (Gelman, p323)
	phi0 <- c(1, alpha.start.values, 1) ## the starting values for the parameters: c(beta0, alpha2-T, gamma2)
	n.params <- length(phi0)

	## sequence of matrices which will converge to DM
	Rt <- array(dim=c(n.params, n.params, max.iter))
	DM <- matrix(nrow=n.params, ncol=n.params)
	DMiter <- rep(0, n.params)

	## sequence of coefficient estimates
	phi.t <- matrix(nrow=n.params, ncol=max.iter+1)
	alpha.t <- alpha.start.values
	phi.t[,1] <- phi0
	loop.idx <- 1:n.params
	while(length(loop.idx)>0) {
		iter <- iter + 1

		##THIS EM CALCULATES PHI^(T+1)
		## E step
		dat <- run.Estep(alpha.t, nlag=nlag,
				 full.data=full.data, assumed.nu=assumed.nu)
		## M step
		phi.t[,iter+1] <- run.Mstep(dat)$phi

		## generate Rt given phi.t
		for(i in loop.idx){
			## defining phi.tmp here as phi^t(i) from Gelman
			phi.t.i <- phi.hat
			phi.t.i[i] <- phi.t[i,iter]
			alpha.t.i <- phi.t.i[2:(n.params-1)]

			## E step
			dat <- run.Estep(alpha.t.i, nlag=nlag,
					 full.data=full.data, assumed.nu=assumed.nu)
			## M step
			phi.tplus1.i <- run.Mstep(dat)$phi

			## fis the ith row of the current Rt matrix
			Rt[i,,iter] <- (phi.tplus1.i-phi.hat)/(phi.t[i,iter]-phi.hat[i])
			if(iter==1) next
			if(all((Rt[i,,iter]-Rt[i,,iter-1])^2 < sqrt(tol))){
				loop.idx <- loop.idx[-which(loop.idx==i)]
				DM[i,] <- Rt[i,,iter]
				DMiter[i] <- iter
			}

		}

		## renew alpha.t
		alpha.tplus1 <- phi.t[2:(n.params-1), iter+1]
		alpha.t <- alpha.tplus1

		if(iter == max.iter) break("maximum iterations reached")

	}

	if(length(loop.idx)>0)	DM <- DM[-loop.idx, -loop.idx]

	DM.out <- list(DM=DM, DMiter=DMiter, loop.idx=loop.idx)
	
	return(DM.out)
}


## E-step of the EM algorithm.

# \arguments{
#  \item{alpha}{Current estimates of the alpha parameters from the GLM model.}
#  \item{full.data}{ A matrix of observed data. See description in EMforCFR helpfile.}
#  \item{nlag}{ The number of time units for lagged data.  Corresponds to length(assumed.nu).}  
#  \item{assumed.nu}{ a vector of probabilities corresponding to the survival distribution, i.e. nu[i]=Pr(surviving i days | fatal case) }
# }

# \value{ A data matrix with the same format as full.data from the EMforCFR() documentation.
# }

run.Estep <- function(alpha, full.data, nlag, assumed.nu){

	## storage for reconstructed deaths
	n.times <- nrow(full.data)/2
	idx.to.reconstruct <- which(full.data[,"grp"]==1 &!is.na(full.data[,"new.times"]))
	times.to.reconstruct <- full.data[idx.to.reconstruct,"time"]
	Dhat.1 <- Dhat.2 <- rep(0, n.times)

	## define the data
	N.1 <- full.data[1:n.times,"N"]
	N.2 <- full.data[(n.times+1):(2*n.times),"N"]
	D.1 <- full.data[1:n.times,"D"]
	D.2 <- full.data[(n.times+1):(2*n.times),"D"]

	## set the alpha coefficients
	alpha.long <- rep(0, n.times)
	alpha.long[times.to.reconstruct[-1]] <- alpha
	alpha.long[(max(times.to.reconstruct)+1):n.times] <- alpha[length(alpha)]

	## ########
	## loops for calculating expected value of the deaths

	for(t in times.to.reconstruct) {
		denom.1 <- denom.2 <- rep(0, nlag)
		for(i in 1:nlag){
			N.idx <- (t+i-1):(t+i-nlag)
			nu.idx <- 1:length(N.idx)
			## group 1
			denom.1[i] <- sum(assumed.nu[nu.idx]*N.1[N.idx]*exp(alpha.long[N.idx]))
			## group 2
			denom.2[i] <- sum(assumed.nu[nu.idx]*N.2[N.idx]*exp(alpha.long[N.idx]))
		}

		## exclude 0s from denominator
		denom.1[denom.1==0] <- NA
		denom.2[denom.2==0] <- NA

		## sum them up
		D.idx <- (t+1):(t+nlag)
		Dhat.1[t] <- N.1[t] * exp(alpha.long[t]) * sum(D.1[D.idx] * assumed.nu[nu.idx]/denom.1, na.rm=TRUE)
		Dhat.2[t] <- N.2[t] * exp(alpha.long[t]) * sum(D.2[D.idx] * assumed.nu[nu.idx]/denom.2, na.rm=TRUE)
	}

	dat <- cbind(Dhat=c(Dhat.1, Dhat.2),
		     new.times=full.data[,"new.times"],
		     grp=full.data[,"grp"],
		     N=full.data[,"N"])
	return(dat)
}


## the M-step

# \arguments{
#  \item{dat}{data matrix passed from EMforCFR().}
# }

# \value{ A list with two components
#         \item{phi }{fitted vector of parameters}
#         \item{Var }{variance-covariance matrix from the fitted model}
# }

run.Mstep <- function(dat){
	subset <- which(dat[,"N"]>0)
	fit <- glm(dat[,"Dhat"] ~ factor(dat[,"new.times"]) + factor(dat[,"grp"]),
		   offset=log(dat[,"N"]),family=poisson, subset=subset)

	phi <- coef(fit)

	out <- list(phi=phi, Var=vcov(fit))
	return(out)
}
