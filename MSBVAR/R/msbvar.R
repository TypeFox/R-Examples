# msbvar() and related functions
# Patrick T. Brandt
# 20081113 : Initial version
# 20120110 : Updated to use new block optimization function by Ryan
#            Davis.  This implements the SWZ block algorithm to get
#            the MLE of the MSVAR and then sets this up for Gibbs
#            sampler.
# 20120112 : Updated to setup some smart starting values and allow the
#            users to input them.  Some error checking is included as
#            well.  Also users multiple random and user input methods
#            to generate starting values for MCMC

# Workhorse msbvar() function with SZ prior to initialize the posterior
# of the Gibbs sampler for the MSBVAR models.

msbvar <- function(Y, z=NULL, p, h,
                   lambda0, lambda1, lambda3,
                   lambda4, lambda5, mu5, mu6, qm,
                   alpha.prior=100*diag(h) + matrix(2, h, h),
                   prior=0, max.iter=40, initialize.opt=NULL)
{
    # Get the number of equations
    m <- ncol(Y)

    # This should be part of a sanity.check.msbvar
    if(h==1)
    {
        stop("\n\n\t -- For MSBVAR models, h>1.  Otherwise, just for a BVAR or VAR!\n")
    }

    # Check the dimensions of alpha.prior
    chk <- dim(alpha.prior)
    if( chk[1]!=h) stop("Incorrect number of rows in alpha.prior.")
    if( chk[2]!=h) stop("Incorrect number of columns in alpha.prior.")

    ########################################################
    # Before the loop this is all initialization for the EM
    # implementation of the Bayesian model.
    ########################################################

    # As default, do a baseline, non-regime model using szbvar() since
    # this gives us all of the inputs we need for later.  User can
    # input their own function / object to do something different as
    # long as it conforms.  See docs for details.

    if(is.null(initialize.opt)==TRUE)  # User provided nothing
    {
        # Inherits all from inputs.
        setup <- initialize.msbvar(Y, p, z, lambda0, lambda1,
                                   lambda3, lambda4, lambda5, mu5,
                                   mu6, nu=m, qm, prior, h,
                                   Q=NULL)

        init.model <- setup$init.model
        Qhat.start <- setup$Qhat.start
        thetahat.start <- setup$thetahat.start
    } else {
    # This is the user input one, so we need to validate it
        if(inherits(initialize.opt$init.model, "BVAR")==FALSE)
        {
            stop("msbvar() initialize.opt list must have an object named init.model of class BVAR.  Create this using szbvar()\n")
        }

        # Check dim of thetahat.start
        tmp <- dim(initialize.opt$thetahat.start)
        if(tmp[1]!=m)
        {
            stop("initialize.opt$thetahat.start has the wrong number of rows\n")
        }
        if(tmp[2]!=(1+m*p+m))
        {
            stop("initialize.opt$thetahat.start has the wrong number of columns\n")
        }
        if(tmp[3]!=h)
        {
            stop("initialize.opt$thetahat.start has the wrong array dimension\n")
        }

       # Validate an initial Q from the user

        if(sum(initialize.opt$Qhat.start)!=h)
        {
            stop("msbvar(): Improper initial Q, transition matrix, for the MS process.  Rows must sum to 1.\n")
        }

        # If we get to here, the input object is valid.  So assign for
        # use
        init.model <- initialize.opt$init.model
        Qhat.start <- initialize.opt$Qhat.start
        thetahat.start <- initialize.opt$thetahat.start

        # Check dim of alpha.prior
        tmp <- dim(alpha.prior)
        if(tmp[1]!=h)
        { stop("msbvar(): Incorrect number of rows in alpha.prior") }

        if(tmp[2]!=h)
        { stop("msbvar(): Incorrect number of columns in alpha.prior")
        }

    }  # End validations of user inputs

    #################################################################
    # EM-blkopt algorithm -- all handled in the Fortran code
    #################################################################

    indms <- 'IAH'   # Right now, user does not have choice over what
                     # switches.  Allow I = Intercepts, A = AR, H =
                     # Variances.

    mlemod <- blkopt(Y, p, thetahat.start, Qhat.start, niter=max.iter,
                     indms)

    # Do setup of the object we want to feed into the posterior
    # sampler, gibbs.msbvar().  This is basically a final pass at
    # setting up the moment matrices for the regressions with the
    # appropriate priors.

    hreg <- setupGibbs(mlemod, Y, init.model)

    # Now re-assign the objects and set up the output

    output <- list(init.model=init.model,
                   hreg=hreg,
                   Q=mlemod$Qhat,
                   fp=mlemod$fpH,
                   m=m, p=p, h=h,
                   alpha.prior=alpha.prior)

    class(output) <- c("MSVARsetup")
    attr(output, "eqnames") <- colnames(Y) # Get variable names for
                                           # attr
    return(output)
}

############################################################
# blkopt() -- the block optimizer for MSBVAR and MS models
############################################################

blkopt <- function(Y, p, thetahat.start, Qhat.start, niter=10, indms) {

    m <- ncol(Y)
    n <- nrow(Y) # later, definition of n will change
    h <- nrow(Qhat.start)

  ##################################################
  ### Setup Regression Matrices
  ##################################################

  # setup LHS Y as Yregmat (Y regression matrix)
    Yregmat <- matrix(Y[(p+1):n,], ncol=m)

  # setup lagged Y's on RHS for AR coefficients, which we define as
  # Xregmat = [Y_t-1, Y_t-2, ..., Y_t-p] (X regression matrix)
  # hence, we will have n-p observations
    Xregmat <- NULL
    if (p>0) { for (i in 1:p) { Xregmat <- cbind(Xregmat, Y[(p-i+1):(n-i),]) } }
    Xregmat <- cbind(rep(1,n-p), Xregmat)  # add intercept on front


  ##################################################
  ### Initialize
  ##################################################

    beta0.it <- array(NA, c(m,1,h))
    beta0.it[,,] <- thetahat.start[,1,]

    betap.it <- NULL
    if (p>0) {
    # rows correspond to equation
    # columns ordered by lag order, then each variable
        betap.it <- array(thetahat.start[,2:(1+p*m),], c(m,m*p,h))
    }

  sig2.it  <- array(thetahat.start[,(1+m*p+1):ncol(thetahat.start),], c(m,m,h))

  theta.it <- array(NA, c(m,1+m*p+m,h))
  theta.it[,1,] <- beta0.it
  if (p>0) theta.it[,2:(1+m*p),] <- betap.it
  theta.it[,(1+m*p+1):ncol(theta.it),] <- sig2.it

  Qhat.it  <- Qhat.start
  # Qhat is in last h columns

  # add one to store llfval for starting values
  llfval <- rep(0,niter+1)

  # run filter to get LLF value for starting values

  # First, obtain residuals
  # e[,,i] is an (n-p) x m matrix containing conditional means
  # for regime i, i=1,2,...,h (e's third dimension is regime)
  e <- array(NA, c(n-p, m, h))
  betahat <- array(NA, c(m,1+m*p,h))
  betahat[,1,] <- beta0.it
  if (p>0) betahat[,2:(1+m*p),] <- betap.it
  # adjust for univariate vs. multivariate case
  if (m>1) {
    for (i in 1:h) { e[,,i] <- Yregmat - Xregmat %*% t(betahat[,,i]) }
  } else {
    for (i in 1:h) { e[,,i] <- Yregmat - Xregmat %*% betahat[,,i] }
  }

  # Using residuals, get filtered regime probabilities
  # HamFilt <- filter.Hamresid(e, sig2.it, Qhat.it)  # R code
  firstHamFilt <- .Fortran("HamiltonFilter",
                           bigt=as.integer(n),
                           m = as.integer(m), p = as.integer(p),
                           h = as.integer(h),
                           e = e,
                           sig2 = sig2.it,
                           Qhat = Qhat.it,
                           f = double(1),
                           filtprSt =
                           matrix(0,as.integer(n-p),as.integer(h))
                           )

  llfval[1] <- firstHamFilt$f

  cat('Initial Log Likelihood Value:', llfval[1] ,'\n', sep=" ")

  ##################################################
  ### Begin Blockwise Optimization
  ##################################################

  cat('Starting blockwise optimization over', niter, 'block optimizations...\n')

  # blocks: intercept(s), AR coefficient(s), variance/Sigma
  for (iter in 1:niter) {
    # intercepts
    beta0.it <- optim(par=c(beta0.it), fn=llf.msar, Y=Yregmat, X=Xregmat, p=p,
                      theta=theta.it, Q=Qhat.it, optstr='beta0',
                      ms.switch=indms, method="BFGS")$par

    if (length(grep('I', indms)) == 0) beta0.it <- array(beta0.it[1:m], c(m,1,h))

    theta.it[,1,] <- array(beta0.it, c(m,1,h))

    # AR coefficients
    if (p>0) {
      betap.it <- optim(par=c(betap.it), fn=llf.msar, Y=Yregmat,
                        X=Xregmat, p=p,
                        theta=theta.it, Q=Qhat.it, optstr='betap',
                        ms.switch=indms, method="BFGS")$par

      # need to check below line...
      if (length(grep('A', indms)) == 0) betap.it <- array(betap.it[1:(m*p)], c(m,m*p,h))
      theta.it[,2:(1+m*p),] <- array(betap.it, c(m,m*p,h))
    }

    # variance/Sigma
    if (m==1) {
      # AR (univariate)
      sig2.it <- optim(par=c(sig2.it), fn=llf.msar, Y=Yregmat, X=Xregmat, p=p,
                       theta=theta.it, Q=Qhat.it, optstr='sig2',
                       ms.switch=indms, method="BFGS")$par

      if (length(grep('H', indms)) == 0) { sig2.it <- array(sig2.it[1], c(m,m,h)) }  # adjust if variance does not switch
      sig2.it <- array(sig2.it, c(m,m,h))
    } else {
      # VAR, only optimize over distinct elements of Sigma for each regime
      sig2.lower <- sig2.it[,,][lower.tri(sig2.it[,,1], diag=TRUE)]
      sig2.lower <- optim(par=c(sig2.lower), fn=llf.msar, Y=Yregmat,
                          X=Xregmat, p=p,
                          theta=theta.it, Q=Qhat.it, optstr='sig2',
                          ms.switch=indms, method="BFGS")$par

      sig2.it  <- array(NA, c(m,m,h))
      # number of distinct: m*(m+1)/2
      nd <- (m*(m+1)/2)
      for (i in 1:h) {
        low <- sig2.lower[(1+(i-1)*nd):(nd+(i-1)*nd)]
        sig2.it[,,i] <- xpnd(low, nrow=m)  # user-defined function below
      }
      if (length(grep('H', indms)) == 0) { sig2.it <- array(sig2.it[,,1], c(m,m,h)) }  # only want first returned value
    }
    theta.it[,(1+m*p+1):ncol(theta.it),] <- sig2.it

    # transition matrix
    Qhat.it <- optim(par=c(Qhat.it[,1:(h-1)]), fn=llf.msar, Y=Yregmat,
                     X=Xregmat, p=p,
                     theta=theta.it, Q=Qhat.it, optstr='Qhat',
                     ms.switch=indms, method="BFGS")$par

    Qhat.it <- matrix(Qhat.it, nrow=h, ncol=h-1)
    Qhat.it <- cbind(Qhat.it, 1-rowSums(Qhat.it))

    # obtain value of likelihood function and store it
    llfval[iter+1] <- -llf.msar(param.opt=beta0.it, Y=Yregmat,
                                X=Xregmat, p=p,
                                theta=theta.it, Q=Qhat.it,
                                optstr='beta0', ms.switch='TOTAL')

    cat('Completed iteration ', iter, ' of ', niter,
        '. Log Likelihood Value: ',
        llfval[iter+1], '\n', sep="")

  }

  # run filter one last time using final parameter estimates

  # First, obtain residuals
  # e[,,i] is an (n-p) x m matrix containing conditional means
  # for regime i, i=1,2,...,h (e's third dimension is regime)
  e <- array(NA, c(n-p, m, h))
  betahat <- array(NA, c(m,1+m*p,h))
  betahat[,1,] <- beta0.it
  if (p>0) betahat[,2:(1+m*p),] <- betap.it
  # adjust for univariate vs. multivariate case
  if (m>1) {
    for (i in 1:h) { e[,,i] <- Yregmat - Xregmat %*% t(betahat[,,i]) }
  } else {
    for (i in 1:h) { e[,,i] <- Yregmat - Xregmat %*% betahat[,,i] }
  }

  # Using residuals, get filtered regime probabilities
  # HamFilt <- filter.Hamresid(e, sig2.it, Qhat.it)  # R code
  lastHamFilt <- .Fortran("HamiltonFilter",
         bigt=as.integer(n),
         m = as.integer(m), p = as.integer(p), h = as.integer(h),
         e = e,
         sig2 = sig2.it,
         Qhat = Qhat.it,
         f = double(1),
         filtprSt = matrix(0,as.integer(n-p),as.integer(h))
         )

  fpH <- lastHamFilt$filtprSt

  output <- list(theta=theta.it,
                 Qhat=Qhat.it,
                 llfval=llfval,
                 fpH=fpH,
                 m=m, p=p, h=h)

  class(output) <- c('blkopt')

  return(output)

}



llf.msar <- function(param.opt, Y, X, p, theta, Q, optstr, ms.switch) {

  m <- ncol(Y)
  n <- nrow(Y) + p
  h <- nrow(Q)

  # initially assign values from theta
  beta0 <- array(theta[,1,], c(m,1,h))
  betap <- NULL
  if (p > 0) betap <- array(theta[,2:(1+m*p),], c(m,m*p,h))
  sig2  <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))
  Qhat  <- Q

  # now choose the parameter over which we are optimizing
  if (optstr=='beta0') {
    beta0 <- array(param.opt, c(m,1,h))
  } else if (optstr=='betap') {
    if (p > 0) betap <- array(param.opt, c(m,m*p,h))
  } else if (optstr=='sig2') {
      sig2  <- array(NA, c(m,m,h))
      # number of distinct: m*(m+1)/2
      nd <- (m*(m+1)/2)
      for (i in 1:h) {
        low <- param.opt[(1+(i-1)*nd):(nd+(i-1)*nd)]
        sig2[,,i] <- xpnd(low, nrow=m)  # user-defined function below
      }
  } else if (optstr=='Qhat') {
    # only passing in first h-1 columns, so add column
    Qhat <- matrix(param.opt, nrow=h, ncol=h-1)
    Qhat <- cbind(Qhat, 1-rowSums(Qhat))
  }

  # numerical checks on Q matrix
  # prevents elements in Q from going negative or greater than 1
  # if that occurs during optimization, then just set to previous Q
  if ( (min(Qhat) <= 0.0001) || (max(Qhat) >= 0.9999 )) Qhat <- Q

  # numerical checks on sig2
  # prevents elements in sig2 from going negative
  # if that occurs during optimization, then just set to previous sig2
  if ( (min(sig2) <= 0.0001) ) sig2 <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))
  # constrain max(off-diagonal) to be less than min(diagonal)
  blnUseOld <- FALSE
  if (m>1) {
    blnSetPast = 0
    for (i in 1:h) {
      if (max(sig2[,,i][lower.tri(sig2[,,i], diag=FALSE)]) > min(diag(sig2[,,i]))) blnSetPast = 1
    }
    if (blnSetPast == 1) blnUseOld <- TRUE
  }
  if (blnUseOld==TRUE) sig2 <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))

  # by default, everything switches, so now adjust by
  # assigning values that do not switch to first state
  # intercept only
  if (ms.switch=='I') {
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='H') { # heteroskedastic
    beta0 <- array(beta0[,,1], c(m,1,h))
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
  } else if (ms.switch=='A') {
    beta0 <- array(beta0[,,1], c(m,1,h))
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='IA') { # homoskedastic
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='IH') {
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
  } else if (ms.switch=='AH') {
    beta0 <- array(beta0[,,1], c(m,1,h))
  }


  #############################################
  # Filtering section
  # Note: there is both Fortran and native R
  #       code in the package (they parallel).
  # Use Fortran for speed, but R for
  # pedagogical purposes.
  #############################################

  # First, obtain residuals
  # e[,,i] is an (n-p) x m matrix containing conditional means
  # for regime i, i=1,2,...,h (e's third dimension is regime)
  e <- array(NA, c(n-p, m, h))
  betahat <- array(NA, c(m,1+m*p,h))
  betahat[,1,] <- beta0
  if (p>0) betahat[,2:(1+m*p),] <- betap
  # adjust for univariate vs. multivariate case
  if (m>1) {
    for (i in 1:h) { e[,,i] <- Y - X %*% t(betahat[,,i]) }
  } else {
    for (i in 1:h) { e[,,i] <- Y - X %*% betahat[,,i] }
  }

  # Using residuals, get filtered regime probabilities
  # HamFilt <- filter.Hamresid(e, sig2.it, Qhat.it)  # R code
  HamFilt <- .Fortran("HamiltonFilter",
         bigt=as.integer(n),
         m = as.integer(m), p = as.integer(p), h = as.integer(h),
         e = e,
         sig2 = sig2,
         Qhat = Qhat,
         f = double(1),
         filtprSt = matrix(0,as.integer(n-p),as.integer(h))
         )

  f <- HamFilt$f

  return(-f) # optim() minimizes negative

}




# This needs to be implemented in C // C++ code later
BHLK.smoother <- function(fp, Q)
{
    TT <- nrow(fp)
    h <- ncol(fp)
    p.smooth <- matrix(0, TT, h)

    p.smooth[TT,] <- fp[TT,]

    for(tt in (TT-1):1)
    {
        p.predict <- Q%*%fp[tt,]
        p.smooth[tt,] <- crossprod(Q, (p.smooth[tt+1,]/p.predict))*fp[tt,]
    }

    return(p.smooth)
}

hregime.setup <- function(h, m, p, TT)
{

    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- matrix(0, 1, h)
    e <- array(0, c(TT, m, h))
    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e))
}


setupGibbs <- function(blkoptobj, Y, init.model)
{

  n <- nrow(Y)
  h <- blkoptobj$h
  m <- blkoptobj$m
  p <- blkoptobj$p

  # now, setup hreg, adjusting for dummies
  hreg <- hregime.reg2.mle(h, m, p, TT=(n-p), fp=blkoptobj$fpH, init.model)


  return(hreg)
}

# Ryan adjusted several items from the original hregime.reg2 function
hregime.reg2.mle <- function(h, m, p, TT, fp, init.model)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- apply(fp, 2, sum)
    e <- array(0, c(TT, m, h))
    Y <- init.model$Y[(m+1+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1+1):nrow(init.model$X),]

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices
        Sxy <- crossprod(X, diag(fp[,i]))%*%Y + crossprod(init.model$X[1:(m+1),], init.model$Y[1:(m+1),])
        Sxx <- crossprod(X, diag(fp[,i]))%*%X + crossprod(init.model$X[1:(m+1),])

        # Compute the regression coefficients
        hstar <- Sxx
        Bk[,,i] <- solve(hstar,Sxy,tol=1e-100)

        # Compute residuals and Sigma (based on Krolzig)

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        Sigmak[,,i] <- (crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]

        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx) #, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}


###########################################################################
# THIS SHOULD BE RENAMED / RE-DEFINED, SINCE EARLIER VERSIONS ARE IN
# THE SVN NOW
###########################################################################
# hregime.reg2
# h  = number of regimes
# m  = number of equations
# p  = number of lags in the VAR
# fp = filter probability matrix, T x h -- where T matches the sample!
# init.model = initial prior / model object produced by szbvar()
###########################################################################

hregime.reg2 <- function(h, m, p, fp, init.model)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- colSums(fp)
    e <- array(0, c(nrow(fp)+m+1, m, h))

    # New version -- just pad in to get the dummy obs.  Use the full
    # design matrix and require first m+1 dummy obs are in each
    # regime.
    Y <- init.model$Y; X <- init.model$X  # So the dummy obs. are
                                          # included here, as defined
                                          # by init.model <- szbvar().

    # Pad the filter probabilities so that each regime gets the dummy
    # observations.

    fp <- rbind(matrix(1, m+1, h), fp)

    # Loop to compute
    # 1) Sums of squares for X and Y for each regime
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices above,
        # so we no longer need to handle these separately as sums
        # inside of each cross-product computation.  See earlier SVN
        # commits to see how this was done (badly) in earlier iterations.

        # Compute the cross products we need..
        # Last one here is so we can avoid ((X'X)^{-1} \ccross Sigma)
        # or the need  to do inverses on potentially sample degenerate
        # values in a regime with too few observations to compute a
        # p.d. version of X'X.

        # Subset the data for each regime
        X1 <- X[fp[,i]==1,]
        Y1 <- Y[fp[,i]==1,]

        Sxy <- crossprod(X1,Y1)
        Sxx <- crossprod(X1)

        ## Sxy <- crossprod(X, diag(fp[,i]))%*%Y
        ## Sxx <- crossprod(X, diag(fp[,i]))%*%X
        ## Syx <- crossprod(Y, diag(fp[,i]))%*%X

        hstar <- Sxx + init.model$H0  # Regression X'X + Prior

        # Compute regressions for a state using QR
        #
        # Note this includes the priors because of the padding in the
        # state-space matrices for fp.

        Bk[,,i] <- qr.coef(qr(hstar), Sxy + init.model$H0[,1:m])

        # Compute residuals and Sigma (based on Krolzig)

        # Get the full residuals -- need these for filtering -- note
        # we do not subset here!
        e[,,i] <- (Y - X%*%Bk[,,i])

        # Now compute the VCOV of the errors.  Note this can be done
        # in several ways for the posterior.  Want to use the
        # regime-specific, sample identified values. This depends on
        # having appropriate df and variation in the error
        # covariance.  Below assumes this is true.  Can do this also
        # with the classic Zellner formula of S0 + Y'Y + H0 - B'(X'X +
        # H0)B, adjusted by df, but this is numerically unstable in
        # some regime classifications.

        # Note this can also handle degenerate regimes (i.e., too few
        # obs.) since the prior will dominate the following
        # computation for the sample variance in a regime and be p.d.

        esub <- as.matrix(e[,,i])
        esub <- esub[fp[,i]==1,]
        Sigmak[,,i] <- (init.model$S0 + init.model$H0[1:m,1:m] +
                        crossprod(esub))/(df[i]+m+1)

        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx)
    }

    # Only return the observation relevant residuals in "e"
    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e[(m+2):nrow(e),,],
                moment=tmp))
}

# Count the number of state transitions for the n(ij) for the
# state-space updates.
count.transitions <- function(s)
  { M <- ncol(s)
    TT <- nrow(s)
    s <- crossprod(t(s), as.matrix(seq(1:M)))
    sw <- matrix(0, M, M)
    for (t in 2:TT)
      { st1 <- s[t-1]
        st <- s[t]
        sw[st1,st] <- sw[st1, st] + 1
      }
    return(sw)
  }
