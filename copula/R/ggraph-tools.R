## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Tools for graphical gof tests based on pairwise Rosenblatt trafo ###########

### FIXME  replace by cCopula() in the long run

##' Conditional copula function C(u2|u1) of u2 given u1
##'
##' @title Bivariate ("2") Conditional Copula function
##' @param u2 data vector (numeric(n))
##' @param u1 data vector (numeric(n))
##' @param family copula family
##' @param theta parameter (for ACs; for elliptical copulas, its rho)
##' @param ... additional args (e.g., df for t copulas)
##' @return C(u2|u1)
##' @author Marius Hofert and Martin Maechler
##' Note: used in Hofert and Maechler (2013)
ccop2 <- function(u2, u1, family, theta, ...) {
    stopifnot(length(u1)==length(u2))
    switch(copFamilyClass(family),
           "ellipCopula"={
               switch(family,
                      "normal"={
                          pnorm((qnorm(u2)-theta*qnorm(u1))/sqrt(1-theta^2))
                      },
                      "t"={
                          stopifnot(hasArg(df))
                          df <- list(...)$df
                          qt.u1 <- qt(u1, df=df)
                          mu <- theta * qt.u1
                          sigma <- sqrt((df+qt.u1^2)*(1-theta^2)/(df+1))
                          pt((qt(u2, df=df)-mu)/sigma, df=df+1)
                      },
                      stop("not yet supported family"))
           },
           "outer_nacopula"={
               ## cacopula(u, cop=onacopulaL(family, list(theta, 1:2)))
               cop <- getAcop(family)
               u <- cbind(u1, u2)
               psiI <- cop@iPsi(u, theta=theta)
	       exp(cop@absdPsi(rowSums(psiI), theta=theta, log=TRUE) -
		   cop@absdPsi(psiI[,1], theta=theta, log=TRUE))
           },
           stop("family ", family, " not yet supported"))
}

##' Compute pairwise Rosenblatt-transformed variables C(u[,i]|u[,j])) for the given
##' matrix u
##'
##' @title Compute pairwise Rosenblatt-transformed variables
##' @param u (n, d)-matrix (typically pseudo-observations, but also perfectly
##'        simulated data)
##' @param cop copula object used for the Rosenblatt transform (H_0;
##'        either outer_nacopula or ellipCopula)
##' @param ... additional arguments passed to ccop2
##' @return (n,d,d)-array cu.u with cu.u[,i,j] containing C(u[,i]|u[,j]) for i!=j
##'         and u[,i] for i=j
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
pairwiseCcop <- function(u, cop, ...)
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2, 0 <= u, u <= 1, d == dim(cop))

    ## 1) determine copula class and compute auxiliary results
    cls <- copClass(cop)
    family <- copFamily(cop) # determine copula family
    switch(cls,
           "ellipCopula"={
               ## build correlation matrix from vector of copula parameters
               P <- diag(1, nrow=d)
	       rho <- cop@parameters
	       P[lower.tri(P)] <- if(family == "normal" || cop@df.fixed)
		   rho else rho[-length(rho)]
               P <- P + t(P)
               diag(P) <- rep.int(1, d)
           },
           "outer_nacopula"={
               ## build "matrix" of dependence parameters
               P <- nacPairthetas(cop)
           },
           stop("not yet supported copula object"))

    ## 2) compute pairwise C(u_i|u_j)
    n <- nrow(u)
    cu.u <- array(NA_real_, dim=c(n,d,d), dimnames=list(C.ui.uj=1:n, ui=1:d, uj=1:d))
    ## cu.u[,i,j] contains C(u[,i]|u[,j]) for i!=j and u[,i] for i=j
    for(i in 1:d) { # first index C(u[,i]|..)
        for(j in 1:d) { # conditioning index C(..|u[,j])
	    cu.u[,i,j] <- if(i==j) u[,i] else
		ccop2(u[,i], u[,j], family, theta=P[i,j], ...)
        }
    }
    cu.u
}


##' Compute matrix of pairwise tests for independence
##'
##' @title Compute matrix of pairwise tests for independence
##' @param cu.u (n,d,d)-array as returned by \code{pairwiseCcop()}
##' @param N number of simulation \code{N} for \code{\link{indepTestSim}()}
##' @param verbose logical indicating if and how much progress info should be printed.
##' @param ... additional arguments passed to indepTestSim()
##' @return (d,d)-matrix of lists with test results (as returned by indepTest())
##'         for the pairwise tests of independence
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
pairwiseIndepTest <-
    function(cu.u, N=256,
	     iTest = indepTestSim(n, p=2, m=2, N=N, verbose = idT.verbose, ...),
	     verbose=TRUE, idT.verbose=verbose, ...)
{
    ## 1) simulate test statistic under independence
    stopifnot(length(dim. <- dim(cu.u)) == 3L)
    stopifnot(dim.[2]==dim.[3])
    n <- dim.[1]
    d <- dim.[2]
    if(verbose)
        cat(sprintf("pairwiseIndepTest( (n=%d, d=%d)): indepTestSim(*, N=%d) .. ",
                    n,d, N))

    ## 2) compute matrix of pairwise p-values for test of independence
    p <- matrix(list(), nrow=d, ncol=d)
    for(j in 1:d) { # column
        if(verbose) if(verbose >= 2) cat("j = ", j, ";  i =", sep="") else cat(j, "")
        uj <- cu.u[,,j] # (n x d)
        for(i in 1:d) { # row
            if(verbose >= 2) cat(i,"")
            p[i,j][[1]] <- if(j==i) list(fisher.pvalue = NA)
	    else indepTest(cbind(uj[,j], uj[,i]), d = iTest)
        }
        if(verbose >= 2) cat("\n")
    }
    if(verbose) cat(" *Sim  done\n")
    p
}

##' Extract p-values from a matrix of indepTest objects
##'
##' @title Extract p-values from a matrix of indepTest objects
##' @param piTest matrix of indepTest objects
##' @return matrix of p-values
##' @author Marius Hofert, Martin Maechler
##' Note: - Kojadinovic, Yan (2010) mention that "fisher.pvalue" performs best;
##'         In d=2 (as we use in pairwiseIndepTest) all three methods are equal.
##'       - used in Hofert and Maechler (2013)
pviTest <- function(piTest){
    matrix(vapply(piTest, function(x) x$fisher.pvalue, numeric(1)),
           ncol=ncol(piTest))
}

##' Computing a global p-value
##'
##' @title Computing a global p-value
##' @param pvalues (matrix of pairwise) p-values
##' @param method vector of methods for p-value adjustments (see ?p.adjust.methods)
##' @param globalFun function determining how to compute a global p-value from a
##'        matrix of pairwise adjusted p-values
##' @return global p-values for each of the specified methods (of how to adjust the
##'         pairwise p-values)
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
gpviTest <- function(pvalues, method=p.adjust.methods, globalFun=min){
    pvalues <- pvalues[!is.na(pvalues)] # vector
    if(all(c("fdr","BH") %in% method))## p.adjust():  "fdr" is alias for "BH"
	method <- method[method != "fdr"]
    sapply(method, function(meth) globalFun(p.adjust(pvalues, method=meth)))
}

## build global pairwise independent test result string
gpviTest0 <- function(pvalues) {
    pvalues <- pvalues[!is.na(pvalues)] # vector
    c("minimum" = min(pvalues),
      "global (Bonferroni/Holm)" = min(p.adjust(pvalues, method="holm")))
}

## string formatter of gviTest0()
gpviString <- function(pvalues, name = "pp-values", sep="   ", digits=2) {
    pv <- gpviTest0(pvalues)
    paste0(name, ":", sep,
           paste(paste(names(pv), sapply(pv, format.pval, digits=digits),
                       sep=": "), collapse=paste0(";",sep)))
}


### Tools for graphical gof tests for copulas with radial parts ################

##' Compute pseudo-observations of the radial part R and the uniform distribution
##' either on the unit sphere (for elliptical copulas) or on the unit simplex (for
##' Archimedean copulas)
##'
##' @title Compute Pseudo-observations of the Radial Part and Uniform Distribution
##' @param x (n, d)-matrix of data
##' @param do.pobs logical indicating whether pseudo-observations should be computed
##' @param method method slot for different copula classes
##' @param ... additional arguments passed to the various methods
##' @return list of two components:
##'         R: n-vector of pseudo-observations of the radial part
##'         S: (n, d)-matrix of pseudo-observations of the uniform distribution on
##'            the unit sphere (for method="ellip") or unit simplex (for method="archm")
##' @author Marius Hofert
##' Note: - the following (true, but unknown) functions can be provided via '...':
##'         "ellip": qQg, the quantile function of the function Q_g in Genest, Hofert, Neslehova (2013)
##'         "archm": iPsi, the inverse of the (assumed) generator
##'       - used in Genest, Hofert, Neslehova (2013)
##'       - for qF for Q-Q plots for "archm", see acR.R
RSpobs <- function(x, do.pobs = TRUE, method = c("ellip", "archm"), ...)
{
    ## check
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    method <- match.arg(method)

    u <- if(do.pobs) pobs(x) else x
    if(!do.pobs && !all(0 <= u, u <= 1))
        stop("'x' must be in the unit hypercube; consider using 'do.pobs=TRUE'")
    switch(method,
           "ellip"={

               ## estimate the standardized dispersion matrix with entries
               ## P_{jk} = \Sigma_{jk}/sqrt{\Sigma_{jj}\Sigma_{kk}}
               ## and compute the inverse of the corresponding Cholesky factor
               ## Note: this is *critical* !!
               ## ----  => completely wrong R's if d > n/2 (roughly)
               P <- as.matrix(nearPD(sin(cor(x, method="kendall")*pi/2),
                                      corr=TRUE)$mat)
	       L <- t(chol(P)) # lower triangular L such that LL' = P
	       ## note: it would be better to stay with 'Matrix' package here and to use LDL

	       ## compute Ys
	       Y <- if(hasArg(qQg)) { # if qQg() has been provided
		   qQg <- list(...)$qQg
		   apply(u, 2, qQg)
	       } else { # estimate via empirical quantiles
                   stop("There is no non-parametric estimator of 'qQg' known yet; provide 'qQg' instead")
                   ## U=(F_1(X_1),..,F_d(X_d)) = (Qg(Y_1),..,Qg(Y_d))
                   ## Our goal is to find the Y's (by applying qQg() to the U's)
                   ## => This is *not* correct:
		   ## x <- scale(x) # standardized data
		   ## sapply(1:ncol(u), function(j)
		   ##        quantile(x[,j], probs=u[,j], names=FALSE))
	       }

	       ## compute Zs and Rs (pseudo-observations of the radial part)
	       ## efficient computation of L^{-1} Y'
	       Z <- solve(L, t(Y))
	       R <- sqrt(colSums(Z^2))

	       ## return R and S (pseudo-obs of uniform distribution on S^d)
	       list(R=R, S=t(Z)/R)

           },
           "archm"={

               if(hasArg(iPsi)) { # if psi^{-1} has been provided
                   iPsi <- list(...)$iPsi
                   iPsiu <- iPsi(u)
                   R <- rowSums(iPsiu)
                   return(list(R=R, S=iPsiu/R)) # return
               }

               ## ... otherwise, estimate estimate psi^{-1} via inverting
               ## the non-parametric estimator of psi of
               ## Genest, Neslehova, Ziegel (2011; Algorithm 1; Section 4.3)

               ## compute r_1,..,r_m and p_1,..,p_m
               wp <- w.p(u) # = w.p(x); p = wp[,"p"]
               d <- ncol(u) # dimension d
               r <- r.solve(wp, d) # compute r_1,..,r_m

               ## compute R_1,..,R_n
               n <- nrow(u)
               R <- R.pobs(r, p=wp[,"p"], n=n) # without shuffling

               ## compute iPsi.n(u) and return
               ## note: we scale the skewed R's (to avoid numerical issues!)
               iPsin <- matrix(iPsi.n(u, r=r/median(R), # scaling
                                      p=wp[,"p"], d=d), ncol=d)
               iPsin.sum <- rowSums(iPsin)
               list(# R = R/med.R, leads to strange discrete shape of Q-Q plots
                    # S = iPsin/(R/med.R), leads to strange discrete shape of Q-Q plots
                    R = iPsin.sum, S = iPsin/iPsin.sum)

           },
           stop("wrong method"))
}

##' @title Compute supposedly Beta distributed observations from supposedly
##'        uniformly distributed observations on the unit sphere
##' @param u (n, d)-matrix of supposedly uniformly distributed (row) vectors
##'        on the unit sphere in IR^d
##' @return (n, d-1)-matrix where the kth column contains supposedly
##'         Beta(k/2, (d-k)/2)-distributed values
##' @author Marius Hofert
##' Note: - see Li, Fang, Zhu (1997); suggestion: take k~=d/2
##'       - used in Genest, Hofert, Neslehova (2013)
gofBTstat <- function(u){
    if (!is.matrix(u)) u <- rbind(u)
    u2 <- u^2
    t(apply(u2, 1, cumsum))[, -ncol(u)]/rowSums(u2)
}

