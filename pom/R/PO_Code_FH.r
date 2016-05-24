siteocc <- function(psi, p, histories, start=NULL, lower=NULL, ...){


require(matrixcalc)

    # -------------  Intermediate functions -------
    get.mod.matrix <- function( f, nsites, nvisits, p.mat ){
    #
    #   Function to return the model matrix from a formula object.  I.e., convert
    #   ~ x1 + x2 into a usable model matrix.
    #
    #   Input:
    #   f = an R formula object, without a response.  I.e., of the form "~ x1 + x2 + ..."
    #   nsites = number of sites = row dimension of matrices
    #   nvisits = number of visits = column dimesion of matrices
    #   p.mat = logical.  If TRUE, assume model in f contains matrices, and return
    #       a 3d array.  Otherwise, assume variables in f are vectors, and return
    #       matrix.
    #
    #   Output:
    #   A matrix with variables in formula as columns.  Also included is
    #   the names of variables in the formula and whether the model has an intercept.
    
        call <- match.call()
        contrasts <- NULL
        mf <- match.call(expand.dots = FALSE)
        mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
        mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
        mf$... <- NULL
        mf$f <- mf$nsites <- mf$nvisits <- mf$p.mat <- NULL

        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")

        mf$formula <- f
	  mf$na.action <- na.pass  # added - FH
        form <- as.character(formula(mf))[2]

        if( nchar(form)==1 & form == "1" ){
            # this is a model with intercept only
            if( p.mat ){
                X <- array(1,c(nsites,nvisits,1))
            } else {
                X <- matrix(1,nsites,1)
            }
            intercept <- TRUE
            nx <- 1
            x.names <- character(0)
            mixture <- "no"
        } else if( form == "Beta.mixture" ){
            #    User is calling for mixture model
            X <- NA
            intercept <- NA
            nx <- 2   # this is number of parameters in beta mixture
            x.names <- character(0)
            mixture <- form
        } else {
            mixture <- "no"
            mf <- eval( mf, parent.frame() )
            mt <- attr( mf, "terms")
            xvars <- as.character(attr(mt, "variables"))[-1]

            if ((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]

            xlev <- if (length(xvars) > 0) {
                    xlev <- lapply(mf[xvars], levels)
                    xlev[!sapply(xlev, is.null)]
            }

            X <- NA
            X <- if ( length(attr( mt, "order")) != 0 ) {
                        model.matrix(mt, mf, contrasts)
                     } else {
                        stop("Empty models not allowed")
                     }

            #   Find out whether intercept is present
            assign.col <- attr( X, "assign" )
            if( sum( assign.col == 0 ) == 1 ){
                intercept <- TRUE
                nx <- sum( unique(assign.col) > 0 )
            } else {
                intercept <- FALSE
                nx <- length( unique(assign.col) )
            }

            x.names <- xvars

            # Add 1 to param count for intercept, if appropriate
            if( intercept ) nx <- nx + 1

            #  Make into a 3-D array if needed
            if( p.mat ){
                if( intercept ){

                    X <- X[,-1]
                    X <- cbind( matrix(1, nsites, nvisits), X)
                }
                dim(X) <- c(nsites, nvisits, nx)
            }

        }

        ans <- list(
            X=X,
            n.covars=nx,
            intercept= intercept,
            vars = x.names,
            mixture=mixture )
        ans
    }

    # -------------------------------------------------

    link <- function(x){
    #   Link function for covariates
        1/(1+exp(-x))   # logistic
    }

    # -------------------------------------------------

     likelihood <- function(beta, X1, X2, hists, k1, k2, P.FUN){
     #  Log likelihood for the patch occupancy problem
     #  Input:
     #  beta = current values of the coefficients, size k1+k2
     #  X1 = 2-D matrix of covariates for psi model, size nsites X k1
     #  X2 = 3-D matrix of covariates for p model, size nsites X nvisits X k2
     #  hists = encounter histories matrix, size nsites X nvisits. values: 0, 1, NA.
     #  k1 = number of coefficients in psi model
     #  k2 = number of coefficients in p model.

          P.FUN <- match.fun( P.FUN )

          psi.coefs = beta[1:k1]
          p.coefs = beta[(k1+1):(k1+k2)]

          # evaluate Psi model
          psi = link(X1%*%psi.coefs)          #psi is nsites X 1

          # evaluate the p part of the likelihood. Sometimes linear, sometimes a mixture
          p.part <- P.FUN( p.coefs, X2, hists)

          k <- as.numeric((rowSums(hists, na.rm=T) >= 1))   # all sites with at least 1 detection
          l <- psi * p.part + (1-psi)*(1-k)

          -sum(log(l))
     }

    # -------------------------------------------------

    beta.mixture <- function( coefs, X, hists){
    #   This implements the beta mixture for the p part of the occupancy likelihood
    #   Note, X must be passed in, but it is not used here.  Just a dummy.
          Ti <- ncol(hists)
          k <- rowSums( hists, na.rm=T )
          b<-coefs[1];
          d<-coefs[2]
          beta(k+b, Ti-k+d)/beta(b, d)
    }

    # -------------------------------------------------

    linear <- function( coefs, X, hists ){
    #   This implments a non-mixture covariate model for the p part of the likelihood.
          nsites <- nrow(hists)
          nvisits <- ncol(hists)

          # Expand p coefficients so multiplication works
          p.coefs <- t(coefs %x% matrix(1, nsites, nvisits))
          dim(p.coefs) <- c(nsites, nvisits, length(coefs))

          p <- X*p.coefs                          # p is nsites X nvisits X k2 here
          p <- apply( p, c(1,2), sum, na.rm=T )   # p is now nsites X nvisits
          p <- link(p)

          p.part <- p^hists * (1-p)^(1-hists)     # p's where hists==NA are wiped out here
          apply( p.part, 1, prod, na.rm=TRUE )
    }


# -------- Main code for siteocc.fn ------
if( missing(histories) ) stop("Site encounter histories matrix must be specified")
if( missing(psi) ) stop("Model for psi (occupancy probabilities) must be specified")
if( missing(p) ) stop("Model for p (encounter probabilties) must be specified")
if(any( unique(unlist(apply(histories,1,unique))) %in% c(0,1,NA) ) == FALSE)stop("Encounter histories must consist of 0's, 1's, and NA's only.")

hist.name <- deparse(substitute(histories))
orig.call <- match.call()
nsites  <- nrow( histories )     #Total # sites sampled
nvisits <- ncol( histories )          #max number of visits to a site

psi.mod.mat <- get.mod.matrix( psi, nsites, nvisits, F )
p.mod.mat   <- get.mod.matrix(   p, nsites, nvisits, T )
X.psi <- psi.mod.mat$X
k.psi <- psi.mod.mat$n.covars   # includes intercept if present

X.p   <- p.mod.mat$X
k.p   <- p.mod.mat$n.covars

#   Decide on mixture and do maximization
if( p.mod.mat$mixture == "Beta.mixture" ){

    # Set up parameters for nlimb function
    if(is.null(start))(start <- rep(0.5, k.psi+k.p))
    if(is.null(lower)){
    	psi.lowlim <- rep(-Inf, k.psi)
    	lower = c(psi.lowlim, 0, 0)
    }
    
    out <- nlminb(start, likelihood, ...,
        X1 = X.psi, X2 = X.p, hists = histories, k1 = k.psi, k2 = k.p, P.FUN = beta.mixture)
    hess <- F.2nd.deriv( out$par, likelihood, X1 = X.psi, X2 = X.p, hists = histories, k1=k.psi, k2=k.p, P.FUN=beta.mixture)  # Compute variance-covariance matrix
    if(is.positive.definite(hess)==FALSE)(cat("Hessian matrix is not positive definite! \nThe Model has not Converged. Parameter Estimates from the last iteration are displayed. \n"))
    vc <- solve(hess)
    se <- sqrt(diag(vc))

    p.cfs <- out$par[(k.psi+1):(k.psi+k.p)]
    se.p <- se[(k.psi+1):(k.psi+k.p)]
    est.p <- choose(ncol(histories),rowSums(histories))*beta.mixture( p.cfs, 1, histories)
    cls <- "mixed.pom"

} else {

    # non-mixture model
    # Set up parameters for nlimb function
    if(is.null(start))(start <- rep(0, k.psi+k.p))
    if(is.null(lower))(lower = -Inf)
    
    out <- nlminb(start, likelihood, ...,
        X1 = X.psi, X2 = X.p, hists = histories, k1 = k.psi, k2 = k.p, P.FUN = linear)
    hess <- F.2nd.deriv( out$par, likelihood, X1 = X.psi, X2 = X.p, hists = histories, k1=k.psi, k2=k.p, P.FUN=linear)  # Compute variance-covariance matrix
    if(is.positive.definite(hess)==FALSE)(cat("Hessian matrix is not positive definite! \nThe Model has not Converged. Parameter Estimates from the last iteration are displayed. \n"))
    vc <- solve(hess)
    se <- sqrt(diag(vc))

    # compute estimated p's
    p.cfs <- out$par[(k.psi+1):(k.psi+k.p)]
    se.p <- se[(k.psi+1):(k.psi+k.p)]
    p.coefs <- t(p.cfs %x% matrix(1, nsites, nvisits))
    dim(p.coefs) <- c(nsites, nvisits, k.p)
    est.p <- X.p*p.coefs                            # p is nsites X nvisits X k2 here
    est.p <- apply( est.p, c(1,2), sum, na.rm=T )   #p is now nsites X nvisits
    est.p <- link(est.p)
    cls <- "pom"
}

#   Compute estimated Psi and p
psi.coefs <- out$par[1:k.psi]
se.psi <- se[1:k.psi]
est.psi <- link(X.psi%*%psi.coefs)          #psi is nsites X 1

#   Find the coefficient names = model terms
#   Assume an intercept is ALWAYS present
psi.names <- as.character(orig.call$psi)[-1]
if( psi.names[1] == "1" ){
	psi.names[1] <- "(Intercept)" 
} else {
	psi.names <- unlist(strsplit(as.character(psi.names), split=" + ", fixed=TRUE))
	psi.names <- unlist(strsplit(as.character(psi.names), split="+", fixed=TRUE))
	psi.names <- c("(Intercept)", psi.names)
}

p.names <- as.character(orig.call$p)[-1]
if( p.names[1] == "1" ){
	p.names[1] <- "(Intercept)" 
} else {
	p.names <- unlist(strsplit(as.character(p.names), split=" + ", fixed=TRUE))
	p.names <- unlist(strsplit(as.character(p.names), split="+", fixed=TRUE))
	p.names <- c("(Intercept)", p.names)
}
if(p.mod.mat$mixture == "Beta.mixture"){
se.p = NULL
p.names <- c("alpha", "beta")
}

names(psi.coefs) <- psi.names
names(p.cfs) <- p.names

ans <- list(    loglik = out$objective,
                convergence = out$convergence,
		convergence.message = out$message,
                call = orig.call,
                naive.psi.est = sum(as.numeric(rowSums(histories, na.rm=T) > 0)) / nsites,
                nsites = nsites,
                nvisits = nvisits,
                psi.coefs = psi.coefs,
                p.coefs = p.cfs,
                se.psi.coefs = se.psi,
                se.p.coefs = se.p,
                hessian = hess,
                psi.ests = est.psi,
                p.ests = est.p,
                aic = 2*out$objective  + 2*(k.psi+k.p),
                bic = 2*out$objective + (k.psi+k.p)*log(nrow(histories)) )
class( ans ) <- cls


ans

}

