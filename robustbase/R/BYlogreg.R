#### http://www.econ.kuleuven.be/public/NDBAE06/programs/roblog/ :
####
#### August 06, 2010  2:14 PM 9121 BYlogreg.r.txt == BYlogreg.R (*this* original)
####    May 04, 2005  9:24 AM 6702 BYlogreg.txt   == BYlogreg.R.~2005~
####    May 04, 2005  9:25 AM 6720 WBYlogreg.txt  == WBYlogreg.R


##  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
##  -------------
##  Christophe Croux, Gentiane Haesbroeck
##  (thanks to Kristel Joossens and Valentin Todorov for improving the code) -
##  ==> Now "contains" both the *weighted* and regular, unweighted  BY-estimator
##
##  This program computes the estimator of Bianco and Yohai in
##  logistic regression. By default, an intercept term is included
##  and p parameters are estimated.
##
##  For more details we refer to
##     Croux, C., and Haesbroeck, G. (2003),
##     ``Implementing the Bianco and Yohai estimator for Logistic Regression'',
##     Computational Statistics and Data Analysis, 44, 273-295
##

## Changes by Martin Maechler, ---> ../man/BYlogreg.Rd
##                                  ------------------

BYlogreg <- function(x0, y, initwml=TRUE, # w.x=NULL,
                     addIntercept=TRUE,
                     const=0.5, kmax = 1000, maxhalf = 10,
                     sigma.min = 1e-4, trace.lev=0)
{
    if(!is.numeric(y))
        y <- as.numeric(y)
    ## if(!is.null(w.x))
    ##     warning("x weights  'w.x'  are not yet made use of")
    if(!is.null(dim(y))) {
        if(ncol(y) != 1)
            stop("y is not onedimensional")
        y <- as.vector(y)
    }
    n <- length(y)

    if(is.data.frame(x0)) {
        x0 <- data.matrix(x0)
    } else if (!is.matrix(x0)) {
        x0 <- matrix(x0, length(x0), 1,
                     dimnames = list(names(x0), deparse(substitute(x0))))
    }
    if(nrow(x0) != n)
        stop("Number of observations in x and y not equal")

    na.x <- !is.finite(rowSums(x0))
    na.y <- !is.finite(y)
    ok <- !(na.x | na.y)
    if(!all(ok)) {
        x0 <- x0[ok, , drop = FALSE]
        y  <- y [ok] # y[ok, , drop = FALSE]
    }

    if(addIntercept) {
        x <- cbind("Intercept" = 1, x0)
    } else { # x0 := x without the  "intercept column"
        x <- x0
        all1 <- apply(x == 1, 2, all)
        if(any(all1))
            x0 <- x[,!all1, drop = FALSE]
        else message("no intercept in the model")
    }

    dx <- dim(x)
    n <- dx[1]
    if(n == 0)
	stop("All observations have missing values!")
    p <- dx[2] # == ncol(x)

    family <- binomial()
    ## Computation of the initial value of the optimization process
    gstart <-
        if(initwml) {
###_ FIXME: Should allow many more schemes:
###_ 1) using MVE with much less singular cases
###_ 2) Instead of {0,1}-weighting with cutoff, w/ weights --> 0 *continuously*
###     --> glm() with "prior" weights instead of 'subset'

            ## hp <- floor(n*(1-0.25))+1
            ##        mcdx <- cov.mcd(x0, quantile.used =hp,method="mcd")
            ##        rdx=sqrt(mahalanobis(x0,center=mcdx$center,cov=mcdx$cov))
            ## mcdx <- CovMcd(x0, alpha=0.75)
            ## rdx  <- sqrt(getDistance(mcdx))
            mcd <- covMcd(x0, alpha=0.75)
            ##                      -----  FIXME: argument!
            D <- sqrt( mahalanobis(mcd$X, mcd$center, mcd$cov) )
            vc  <- sqrt(qchisq(0.975, p-1))
            ##                 -----       FIXME: 'vc' should be argument!
            wrd <- D <= vc
### FIXME_2: use  weights and "weights.on.x'  as in Mqle ( ./glmrobMqle.R )
	    ## glm(y~x0, family=binomial, subset = wrd)$coef
	    glm.fit(x[wrd,,drop=FALSE], y[wrd], family=family)$coef
	} else {
	    glm.fit(x, y, family=family)$coef
        }

    sigma1 <- 1/sqrt(sum(gstart^2))
    xistart <- gstart*sigma1
    stscores <- x %*% xistart

    ## Initial value for the objective function
    oldobj <- mean(phiBY3(stscores/sigma1, y, const))

    converged <- FALSE
    kstep <- 1L
    while(kstep < kmax && !converged)
    {
        unisig <- function(sigma) mean(phiBY3(stscores/sigma, y, const))
        ##                             ------
        optimsig <- nlminb(sigma1, unisig, lower=0)# "FIXME" arguments to nlminb()
        ##          ======
        if(trace.lev) cat(sprintf("k=%2d, s1=%12.8g: => new s1= %12.8g",
                                  kstep, sigma1, optimsig$par))# MM: jhalf =!?= 1 here ??
        sigma1 <- optimsig$par

        if(sigma1 < sigma.min) {
            if(trace.lev) cat("\n")
            warning(gettextf("Implosion: sigma1=%g became too small", sigma1))
            kstep <- kmax #-> *no* convergence
        } else {
            ## gamma1 <- xistart/sigma1
            scores <- stscores/sigma1
            newobj <- mean(phiBY3(scores, y,const))
            oldobj <- newobj
            grad.BY <- colMeans((derphiBY3(scores,y,const) %*% matrix(1,ncol=p))*x)
            h <- -grad.BY + (grad.BY %*% xistart) *xistart
            finalstep <- h/sqrt(sum(h^2))

            if(trace.lev) {
                if(trace.lev >= 2) cat(sprintf(", obj=%12.9g: ", oldobj))
                cat("\n")
            }

            ## FIXME repeat { ... }   {{next 4 lines are also inside while(..) below}}
            xi1 <- xistart+finalstep
            xi1 <- xi1/sum(xi1^2)
            scores1 <- (x %*% xi1)/sigma1
            newobj <- mean(phiBY3(scores1,y,const))

            ## If 'newobj' is not better, try taking a smaller step size:
            hstep <- 1.
            jhalf <- 1L
            while(jhalf <= maxhalf & newobj > oldobj)
            {
                hstep <- hstep/2
                xi1 <- xistart+finalstep*hstep
                xi1 <- xi1/sqrt(sum(xi1^2))
                scores1 <- x %*% xi1/sigma1
                newobj <- mean(phiBY3(scores1,y,const))
                if(trace.lev >= 2)
                    cat(sprintf("  jh=%2d, hstep=%13.8g => new obj=%13.9g\n",
                                jhalf, hstep, newobj))
                jhalf <- jhalf+1L
            }

            converged <-
                not.improved <- (jhalf > maxhalf && newobj > oldobj)
            if(not.improved) {
                ## newobj is "worse" and step halving did not improve
                message("Convergence Achieved")
            } else {
                jhalf <- 1L
                xistart <- xi1
                oldobj <- newobj
                stscores <- x %*% xi1
                kstep <- kstep+1L
            }
        }
    } ## while( kstep )

    if(kstep == kmax) {
        warning("No convergence in ", kstep, " steps.")
        list(convergence=FALSE, objective=0, coefficients= rep(NA,p))
    } else {
        gammaest <- xistart/sigma1
        V <- vcovBY3(x, y, const, estim=gammaest, addIntercept=FALSE)
        list(convergence=TRUE, objective=oldobj, coefficients=gammaest,
             cov = V, sterror = sqrt(diag(V)),
             iter = kstep)
    }
}


### -- FIXME: nlminb() allows many tweaks !!
### -- -----  but we use nlminb() for ONE-dim. minimization over { sigma >= 0 } - really??
##  MM: my version would rather use  optimize() over over  log(sigma)
glmrobBY.control <-
    function(maxit = 1000, const = 0.5, maxhalf = 10)
    ## FIXME: sigma.min
    ## MM: 'acc' seems a misnomer to me, but it's inherited from  MASS::rlm
    ## TODO acc = 1e-04, test.acc = "coef", tcc = 1.345)
{
    ## if (!is.numeric(acc) || acc <= 0)
    ##     stop("value of acc must be > 0")
    ## if (test.acc != "coef")
    ##     stop("Only 'test.acc = \"coef\"' is currently implemented")
    ## if (!(any(test.vec == c("coef", "resid"))))
    ##	  stop("invalid argument for test.acc")
    if(!is.numeric(maxit) || maxit <= 0)
	stop("maximum number of \"kstep\" iterations must be > 0")
    if(!is.numeric(maxhalf) || maxhalf <= 0)
	stop("maximal number of *inner* step halvings must be > 0")
    ## if (!is.numeric(tcc) || tcc <= 0)
    ##     stop("value of the tuning constant c (tcc) must be > 0")
    if(!is.numeric(const) || const <= 0)
	stop("value of the tuning constant c ('const') must be > 0")

   list(## acc = acc, consttest.acc = test.acc,
         const=const,
         maxhalf=maxhalf,
         maxit=maxit #, tcc = tcc
         )
}


##' @param intercept logical, if true, X[,] has an intercept column which should
##'                  not be used for rob.wts
glmrobBY <- function(X, y,
                     weights = NULL, start = NULL, offset = NULL,
                     method = c("WBY","BY"), weights.on.x = "none",
                     control = glmrobBY.control(...), intercept = TRUE,
                     trace.lev = 0, ...)
{
### THIS is *NOT* exported

    method <- match.arg(method)
    if(!is.null(weights) || any(weights != 1)) ## FIXME (?)
        stop("non-trivial prior 'weights' are not yet implemented for \"BY\"")
    if(!is.null(start))
        stop(" 'start' cannot yet be passed to glmrobBY()")
    if(!is.null(offset))
        stop(" 'offset' is not yet implemented for \"BY\"")
    const   <- if(is.null(cc <- control$const  )) 0.5 else cc
    kmax    <- if(is.null(cc <- control$maxit  )) 1e3 else cc
    maxhalf <- if(is.null(cc <- control$maxhalf))  10 else cc
    if(!identical(weights.on.x, "none"))
        stop("'weights.on.x' = ", format(weights.on.x)," is not implemented")
    ## w.x <- robXweights(weights.on.x, X=X, intercept=intercept)
    ##
    ## MM: all(?) the  BY3() functions below would need to work with weights...

    r <- BYlogreg(x0=X, y=y, initwml = (method == "WBY"), ## w.x=w.x,
		  addIntercept = !intercept, ## add intercept if there is none
		  const=const, kmax=kmax, maxhalf=maxhalf,
		  ## FIXME sigma.min  (is currently x-scale dependent !????)
		  trace.lev=trace.lev)
    ## FIXME: make result more "compatible" with other glmrob() methods
    r
}



### Functions needed for the computation of estimator of Bianco and Yohai ----------------------

## From their paper:

## A last remark is worth mentioning: when huge outliers occur in
## the logistic regression setting, often numerical imprecision occurs in the computation
## of the deviances given by
##    d(s;y_i)= -y_i log F(s) - (1-y_i) log{1-F(s)} .
##
## Instead of directly computing this expression, it can be seen that a
## numerically more stable and accurate formula is given by
##    log(1 + exp(-abs(s))) + abs(s)* ((y-0.5)*s < 0)
## in which the second term equals abs(s) if the observation is misclassified, 0 otherwise.
dev1 <- function(s,y) log(1+exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev2 <- function(s,y) log1p(exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev3 <- function(s,y) -( y  * plogis(s, log.p=TRUE) +
                        (1-y)*plogis(s, lower.tail=FALSE, log.p=TRUE))
## MM[FIXME]: first tests indicate that  dev3() is clearly more accurate than
##            their dev1() !!
## MM{FIXME2}: In code below have (or "had") three cases of same formula, but
##             with 's>0' instead of 's<0' : This is == dev?(-s, y) !!

## for now,  100% back-compatibility:
devBY <- dev1
rm(dev1, dev2, dev3)

## MM: This is from my vignette, but  *not* used
log1pexp <- function(x) {
    if(has.na <- any(ina <- is.na(x))) {
	y <- x
	x <- x[ok <- !ina]
    }
    t1 <- x <= 18
    t2 <- !t1 & (tt <- x <= 33.3)
    r <- x
    r[ t1] <- log1p(exp(x[t1]))
    r[ t2] <- { x2 <- x[t2]; x2 + exp(-x2) }
    r[!tt] <- x[!tt]
    if(has.na) { y[ok] <- r ; y } else r
}


phiBY3 <- function(s,y,c3)
{
  s <- as.double(s)
  ## MM FIXME  log(1 + exp(-.)) ... but read the note above !! ---
  dev. <- devBY(s,y)
  ## FIXME: GBY3Fs()  computes the 'dev' above *again*, and
  ##        GBY3Fsm() does with 's>0' instead of 's<0'
  rhoBY3(dev.,c3) + GBY3Fs(s,c3) + GBY3Fsm(s,c3)
}

rhoBY3 <- function(t,c3)
{
    ec3 <- exp(-sqrt(c3))
    t*ec3* (t <= c3) +
        (ec3*(2+(2*sqrt(c3))+c3) - 2*exp(-sqrt(t))*(1+sqrt(t)))* (t > c3)
}

psiBY3 <- function(t,c3)
{
    exp(-sqrt(c3)) *(t <= c3) +
    exp(-sqrt( t)) *(t >  c3)
}
## MM: This is shorter (but possibly slower when most t are <= c3 :
## psiBY3 <- function(t,c3) exp(-sqrt(pmax(t, c3)))

##'   d/dt psi(t, c3)
derpsiBY3 <- function(t, c3)
{
    r <- t
    r[in. <- (t <= c3)] <- 0
    if(any(out <- !in.)) {
        t <- t[out]
        st <- sqrt(t)
        r[out] <- -exp(-st)/(2*st)
    }
    r
}

## MM: FIXME  this is not used above
sigmaBY3 <- function(sigma,s,y,c3)
{
    mean(phiBY3(s/sigma,y,c3))
}

derphiBY3 <- function(s,y,c3)
{
    Fs <- exp(-devBY(s,1))
    ds <- Fs*(1-Fs) ## MM FIXME: use expm1()
    dev. <- devBY(s,y)
    Gprim1 <- devBY(s,1)
    Gprim2 <- devBY(-s,1)

    -psiBY3(dev.,c3)*(y-Fs) + ds*(psiBY3(Gprim1,c3) - psiBY3(Gprim2,c3))
}

der2phiBY3 <- function(s, y, c3)
{
    s <- as.double(s)
    Fs <- exp(-devBY(s,1))
    ds <- Fs*(1-Fs) ## MM FIXME: use expm1()
    dev. <- devBY(s,y)
    Gprim1 <- devBY(s,1)
    Gprim2 <- devBY(-s,1)
    der2 <- derpsiBY3(dev.,c3)*(Fs-y)^2  + ds*psiBY3(dev.,c3)
    der2 <- der2+ ds*(1-2*Fs)*(psiBY3(Gprim1,c3) - psiBY3(Gprim2,c3))
    der2 - ds*(derpsiBY3(Gprim1,c3)*(1-Fs) +
               derpsiBY3(Gprim2,c3)*  Fs )
}


GBY3Fs <- function(s,c3)
{
    e.f <- exp(0.25)*sqrt(pi)
    ## MM FIXME: Fs = exp(..) and below use  log(Fs) !!
    Fs <- exp(-devBY(s,1))
    resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
     ## MM FIXME: use expm1():
    resGinf <- (resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
    resGsup <- ((Fs*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
        as.numeric(s > -log(exp(c3)-1))
    resGinf + resGsup
}


GBY3Fsm <- function(s,c3)
{
    e.f <- exp(0.25)*sqrt(pi)
    ## MM FIXME: Fsm = exp(..) and below use  log(Fsm) !!
    Fsm <- exp(-devBY(-s,1))
    resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
     ## MM FIXME: use expm1():
    resGinf <- (resGinf+(Fsm*exp(-sqrt(-log(Fsm))))) * as.numeric(s >= log(exp(c3)-1))
    resGsup <- ((Fsm*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
        as.numeric(s < log(exp(c3)-1))
    resGinf + resGsup
}

## Compute the standard erros of the estimates -
##  this is done by estimating the asymptotic variance of the normal
##  limiting distribution of the BY estimator - as derived in Bianco
##  and Yohai (1996)
##
sterby3 <- function(x0, y, const, estim, addIntercept)
{
    sqrt(diag(vcovBY3(x0, y, const=const, estim=estim, addIntercept=addIntercept)))
}

vcovBY3 <- function(z, y, const, estim, addIntercept)
{
    stopifnot(length(dim(z)) == 2)
    if(addIntercept) z <- cbind(1, z)
    d <- dim(z)
    n <- d[1]
    p <- d[2]
    argum <- z %*% estim
    matM <- IFsqr <- matrix(0, p, p)
    for(i in 1:n)
    {
        myscalar <- as.numeric(der2phiBY3(argum[i],y[i], c3=const))
        zzt <- tcrossprod(z[i,])
        matM <- matM + myscalar * zzt
        IFsqr <- IFsqr + derphiBY3(argum[i],y[i], c3=const)^2 * zzt
    }

    matM    <- matM/n
    matMinv <- solve(matM)
    IFsqr <- IFsqr/n
    ## Now,  asymp.cov  =  matMinv %*% IFsqr %*% t(matMinv)

    ## provide  vcov(): the full matrix
    (matMinv %*% IFsqr %*% t(matMinv))/n
}


