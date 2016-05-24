## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### is source()d from  demo(estimation.gof)  and also in some of the tests
##		       ------------------- ==> ../../demo/estimation.gof.R
##					       ~~~~~~~~~~~~~~~~~~~~~~~~~~~

##' measures user run time in milliseconds
utms <- function(x) 1000 * system.time(x)[[1]]
##' formats such utms() times "uniformly":
f.tms <- function(x) paste(round(x),"ms") # with a space (sep=" ") !

##' Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
##' copulas
##'
##' @title Fitting and Goodness-Of-Fit for Archimedean copulas
##' @param n sample size
##' @param d dimension
##' @param simFamily *Archimedean* family to be sampled
##' @param tau degree of dependence of the sampled family in terms of Kendall's tau
##' @param N number of bootstrap replications
##' @param include.K whether or not K is included
##' @param n.MC if > 0 it denotes the sample size for SMLE
##' @param estim.method  estimation method (see fitCopula or gofCopula)
##' @param gof.method goodness-of-fit transformation (see gofCopula)
##' @param checkFamilies vector of Archimedean families to be used for gof
##' @param verbose
##' @return a numeric matrix ...
##' @author Marius Hofert and Martin Maechler
##' @note In a next step, one could include Ga and t as well, not only Archimedean copulas
estimation.gof <- function(n, d, simFamily, tau,
			   N = 1, # dummy number of bootstrap replications
                           estim.method = eval(formals(fitCopula)$method),
                           simulation="pb", verbose=TRUE,
                           gof.trafo = c("rtrafo", "htrafo"),
                           gof.method = eval(formals(gofTstat)$method),
			   checkFamilies = setdiff(.ac.longNames, "AMH")) # as AMH is not available for d > 2
{
    ## generate data
    theta <- getAcop(simFamily)@iTau(tau)
    simCop <- onacopulaL(simFamily, list(theta, 1:d)) # not possible with archmCopula() as AMH is not availble for d > 2
    if(verbose){
        r <- function(x) round(x,4) # for output
	cat("\n\n### Output for estim.method = \"",estim.method, "\" gof.trafo = \"",
            gof.trafo,"\", and gof.method = \"",gof.method,"\"\n\n",sep="")
    }
    U <- rCopula(n, simCop)

    ## estimation and gof
    estim.method <- match.arg(estim.method)
    gof.trafo <- match.arg(gof.trafo)
    gof.method <- match.arg(gof.method)
    n.cf <- length(checkFamilies)
    sig <- logical(n.cf)
    tau <- gof <- ute <- utg <- est <- numeric(n.cf)
    for(k in seq_len(n.cf)) {
        ## estimate the parameter with *fitCopula* and the provided method
        cop <- archmCopula(family=checkFamilies[k], dim=d) # use archmCopula() due to fitCopula()
        if(verbose)
            cat("Estimation and GoF for ", checkFamilies[k], ":\n\n", sep="")
	ute[k] <- utms(fit <- fitCopula(cop, data=U, method=estim.method,
                                        estimate.variance=FALSE))
        est[k] <- fit@estimate
        tau[k] <- tau(fit@copula)
        if(verbose){
            cat("   theta hat      = ",r(est[k]),"\n",
                "   tau hat        = ",r(tau[k]),"\n",
                ## The exact string 'Time ' must appear at the beginning of line
                ## for 'R CMD Rdiff' to *not* look at differences there:
                "Time estimation   = ",f.tms(ute[k]),"\n", sep="")
	}
        cop@parameters <- est[k] # set the parameter of the estimated copula; or use: fit@copula
        ## apply a goodness-of-fit test to the estimated copula
        ## {{ use rtrafo() or htrafo() if you want the transformed u }}
	utg[k] <-
	    utms(gof[k] <-
		 gofCopula(cop, x=U, N=N, method=gof.method,
			   estim.method=estim.method, simulation=simulation,
                           verbose=verbose,
                           trafo.method=gof.trafo)$p.value)
	sig[k] <- (gof[k] < 0.05) # TRUE/FALSE <--> 1/0 -- Careful: may be NA !
	if(verbose)
	    cat("   p-value	   = ",r(gof[k]),"\n",
		"   < 0.05	   = ", format(sig[k]), "\n",
		"Time GoF comp = ",f.tms(utg[k]),"\n\n", sep="")
    }

    ## result
    names(est) <- checkFamilies
    cbind(theta_hat=est, tau_hat  =tau,
	  timeEstim=ute,
	  P_value  =gof, "< 0.05" =sig,
	  timeGoF  =utg)
}
