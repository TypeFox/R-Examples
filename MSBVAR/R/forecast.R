# Generate forecast methods for the models in the MSBVAR package
#
# Patrick T. Brandt
#
# 20120621 : Updated to include MSBVAR forecast functions.

"forecast" <- function(varobj, nsteps, A0=t(chol(varobj$mean.S)),
                       shocks=matrix(0,nrow=nsteps,ncol=dim(varobj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps,ncol=nrow(varobj$exog.coefs)),
                       N1, N2)
{
    if(inherits(varobj,"VAR")){
        return(forecast.VAR(varobj, nsteps, A0=A0,
                            shocks=shocks, exog.fut=exog.fut))
    }

    if(inherits(varobj, "BVAR")){
        return(forecast.VAR(varobj, nsteps, A0=A0,
                            shocks=shocks, exog.fut=exog.fut))
    }

    if(inherits(varobj, "BSVAR")){
        return(forecast.VAR(varobj, nsteps, A0=solve(varobj$A0.mode),
                       shocks=shocks, exog.fut=exog.fut))
    }

    if(inherits(varobj, "MSBVAR")){
        return(forecast.MSBVAR(x=varobj, k=nsteps, N1, N2))
    }
}

# This is the generic VAR forecasting function.  The other
"forecast.VAR" <-
function(varobj, nsteps, A0, shocks, exog.fut)
{
   # Set up the initial parameters for the VAR forecast function from
   #  VAR object
  y <- varobj$y
  intercept <- varobj$intercept
  ar.coefs <- varobj$ar.coefs
  exog.coefs <- varobj$exog.coefs
  m<-dim(ar.coefs)[1]
  p<-dim(ar.coefs)[3]
  capT<-nrow(y)
  yhat<-rbind(y,matrix(0,ncol=m,nrow=nsteps))

   # Compute the deterministic part of the forecasts (less the intercept!)
   if(is.na(sum(varobj$exog.coefs))==F)
     {
       deterministic.VAR <- as.matrix(exog.fut) %*% exog.coefs
     }
   else
     { deterministic.VAR <- matrix(0,nrow=nsteps,ncol=m)
     }

   # Now loop over the forecast horizon
   for(h in 1:nsteps)
     {  yhat[capT + h, ] <- (yhat[capT + h - 1,] %*% ar.coefs[,,1] +
                             intercept + deterministic.VAR[h,] + (shocks[h,]%*%A0))
       if (p>1) {for(i in 2:p)
       { yhat[capT + h, ] <- (yhat[capT + h, ] +
                              (yhat[capT + h - i, ] %*% ar.coefs[,,i]))

       }}
     }
  output <- ts(yhat, start = start(varobj$y), frequency = frequency(varobj$y), names = colnames(varobj$y))
  attr(output, "class") <- c("forecast.VAR", "mts", "ts")
  attr(output, "eqnames") <- attr(varobj, "eqnames")
  return(output)
}

"forecast.BVAR" <- function(varobj, nsteps, A0, shocks, exog.fut)
{
    output <- forecast.VAR(varobj, nsteps, A0, shocks, exog.fut)
    attr(output, "class") <- c("forecast.BVAR", "mts", "ts")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"forecast.BSVAR" <- function(varobj, nsteps, A0=solve(varobj$A0.mode), shocks, exog.fut)
{
    output <- forecast.VAR(varobj, nsteps, A0, shocks, exog.fut)
    attr(output, "class") <- c("forecast.BSVAR", "mts", "ts")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"uc.forecast" <- function(varobj, nsteps, burnin, gibbs,
                          exog=NULL)
{
    if(inherits(varobj, "VAR"))
    {
        stop("Not implemented for VAR models!\nUse a BVAR with a flat-flat prior if you want this case.\n")
##         varobj$H0 <- matrix(0, nrow(varobj$Bhat), nrow(varobj$Bhat))
##         varobj$S0 <- matrix(0, ncol(varobj$Bhat), ncol(varobj$Bhat))
##         output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs, exog)
##         attr(output, "class") <- c("forecast.VAR")
##         return(output)
    }

    if(inherits(varobj, "BVAR"))
    {
        output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BSVAR"))
    {
      stop("Not yet implemented for BSVAR models!\n")
##       output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs,exog)
##       attr(output, "class") <- c("uc.forecast.VAR", "mts", "ts")
##       return(output)
    }
}

"uc.forecast.VAR" <- function(varobj, nsteps, burnin, gibbs, exog)
  { # Extract all the elements from the VAR object
    y <- varobj$y
    ar.coefs <- varobj$ar.coefs
    intercept <- varobj$intercept
    A0 <- t(chol(varobj$mean.S))
    X <- varobj$X         # rhs variables for the model
    Y <- varobj$Y         # lhs variables for the model
    H0 <- varobj$H0
    S0 <- varobj$S0
#    mu <- varobj$hyperp
    exog.coefs <- varobj$exog.coefs
    z <- varobj$z
#    lambda0 <- varobj$prior[1]
#    lambda1 <- varobj$prior[2]
#    lambda3 <- varobj$prior[3]
#    lambda4 <- varobj$prior[4]
#    lambda5 <- varobj$prior[5]
#    mu5 <- varobj$prior[6]
#    mu6 <- varobj$prior[7]
    nu <- varobj$prior[8]
    prior <- varobj$prior.type
    num.exog <- varobj$num.exog
    qm <- varobj$qm
    ncoef <- nrow(varobj$Bhat)

  # Get some constants we are going to need
    starttime<-date()           # Starting time for simulation

    p<-dim(ar.coefs)[3]         # Capture the number of lags
                                # from input ar coefficients

    m<-ncol(y);                 # Number of endogenous
                                # variables in the VAR
    k<-m*nsteps;              # k = mh, the maximal number
                              # of forecasts

    capT<-nrow(y)               # Number of observations we
                                # are going to use.

  # Make arrays to hold the Gibbs sampler results
  yforc<-array(0,c(gibbs,nsteps,m))

  # Do the Gibbs draws.....
  for(i in 1:(burnin+gibbs))
    { # Step (a): Compute a draw of the conditional forecasts
      # COMPUTE INNOVATIONS: These are the structural innovations

      # First draw the innovations
      epsilon.i <- matrix(rnorm(nsteps*m),nrow=nsteps,ncol=m)

      # Then construct a forecast using the innovations
      ytmp <- forecast.VAR(varobj, nsteps, A0=A0, shocks=epsilon.i,
                           exog)

      # Store draws that are past the burnin in the array
      if(i>burnin)
        { for(j in 1:m)
            { yforc[(i-burnin),(1:nsteps),j]<-ytmp[((capT+1):(capT+nsteps)),j] }
        }


      # Step (b): Compute the mode of the posterior for the
      # forecast distribution.  This is the "extended"
      # dataset that includes the i'th Gibbs sample forecast

      # Set up the updated Y Matrix
      # this is just "ytmp" from above

      Y.update <- ytmp[(capT-p+1):nrow(ytmp),]

      # Set up the updated X -- this is hard because we need to get
      # the RHS lags correct.  We do this by padding the existing Y
      # and then building the lags.  This reuses the code for the lag
      # construction in the szbvar code

      # Now build rhs -- start with an empty matrix
      X.update <- matrix(0, nsteps, m*p+1)
        # Add in the constant
      X.update[, (m*p+1)] <- matrix(1, nsteps, 1)
        # Put in the lagged y's

      # Note that constant is put last here when the lags are made
      for(j in 1:p)
        {
          X.update[1:nsteps,(m*(j-1)+1):(m*j)] <- matrix(Y.update[(p+1-j):(nsteps+p-j),],ncol=m)
        }

      # Put in exogenous coefficients if there are any.
      if(is.null(exog)==F)
        {
          X.update<-cbind(X.update,exog);
        }

      # Now, stack the original Y data and the augmented data.
      Y.update <- rbind(Y, ytmp[(capT+1):nrow(ytmp),])

      # Set up crossproducts and inverses we need
      X.update <- rbind(X, X.update)
      XX.update <- crossprod(X.update)     # Cross product of RHS variables
      hstar.update <- H0 + XX.update       # Prior + Cross product

      # Updated Regression estimates, Beta | Sigma
      B.update<-solve((hstar.update),(crossprod(X.update,Y.update) + H0[,1:m]))

      # Posterior mean of Sigma | Beta
      S.update <- (S0 + crossprod(Y.update)
                   + H0[1:m,1:m] - t(B.update)%*%(hstar.update)%*%(B.update))/(capT+nsteps+nu-m-1)

      # Posterior variance of B: VBh = diag(Sh.*.(inv((H0 + x'x)))
      hstarinv <- solve(hstar.update)
      vcv.Bh <- kronecker(S.update,hstarinv)

      # Draw from the conditional posterior pdfs of the parameters

      # This is only valid for just-identified models.
      df <- capT - m*p - m - 1 + nsteps
      wisharts <- rwishart(1, df, diag(m))

      # Generate the draws from the Wishart and the Beta
      # Wishart draw
      Sigmat <- (chol(S.update))
      Sigma.Draw <- t(Sigmat)%*%(df*solve(matrix(wisharts,m,m)))%*%Sigmat
      sqrtwish <- t(chol(Sigma.Draw))
      # Covariance of beta
      bcoefs.covar <- t(chol(vcv.Bh))

      # Draw beta|Sigma
      aplus <- matrix(B.update, ncol=1) + bcoefs.covar%*%matrix(rnorm(nrow(bcoefs.covar)), ncol=1)
      aplus <- matrix(aplus, ncol=m)

      aplus.coefs<-t(aplus[1:(m*p),]);        # extract the ar coefficients
      dim(aplus.coefs)<-c(m,m,p)                    # push ar coefs into M x M x P array
      aplus.coefs<-aperm(aplus.coefs,c(2,1,3))      # reorder array


      intercept <- aplus[(m*p+1),]       # get drawn intercept....
      ar.coefs<-aplus.coefs            # AR coefs
      A0 <- sqrtwish

      if(num.exog!=0)
        {
          exog.coefs <- aplus[(m*p+2):nrow(aplus),]
        }


      # Print some intermediate results to capture progress....
      # and tell us that things are still running
      if (i%%500==0)
        { cat("Gibbs Iteration = ", i, "     \n");
          if(i<=burnin)
            { cat("(Still a burn-in draw.)\n");
            }

        }
      # Back to the top of the Gibbs loop....
    }
  endtime<-date()
  # Print time stamp so we know how long everything took.
  cat("Start time : ", starttime, "\n");
  cat("End time   : ", endtime, "\n");
    # Returns a list object
    output <- list(forecast=yforc)
#    attr(output, "class") <- c("forecast.VAR")
    return(output)
}

"hc.forecast" <- function(varobj, yconst, nsteps, burnin, gibbs, exog=NULL)
{
    if(inherits(varobj, "VAR"))
    {
        stop("Not yet implemented for VAR models!\nUse a BVAR model.")
##         output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
##                                   gibbs, exog)
##         attr(output, "class") <- c("forecast.VAR")
##         return(output)
    }

    if(inherits(varobj, "BVAR"))
    {
        output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
                                  gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BSVAR"))
    {
      stop("Not yet implemented for B-SVAR models!\n")
##         output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
##                                   gibbs, exog)
##         attr(output, "class") <- c("hc.forecast.VAR", "mts", "ts")
##         return(output)
    }
}

"hc.forecast.VAR" <-
function(varobj, yconst, nsteps, burnin, gibbs, exog=NULL)
{
    # Extract all the elements from the VAR object that we will need
    y <- varobj$y
    ar.coefs <- varobj$ar.coefs
    intercept <- varobj$intercept
    exog.coefs <- varobj$exog.coefs
    A0 <- t(chol(varobj$mean.S))
    mu <- varobj$hyperp
    prior <- varobj$prior
    ncoef <- nrow(varobj$Bhat)
    X <- varobj$X         # rhs variables for the model
    Y <- varobj$Y         # lhs variables for the model
    H0 <- varobj$H0       # precision for the var coefs
    S0 <- varobj$S0       # precision for the Sigma
    nu <- varobj$prior[8] # df

    # Get some constants we are going to need from the inputs
    starttime<-date()           # Starting time for simulation
                                # of forecasts

    q<-nrow(as.matrix(yconst));   # Number of restrictions

    p<-dim(ar.coefs)[3]         # Capture the number of lags
                                # from input ar coefficients

    m<-ncol(y);                 # Number of endogenous
                                # variables in the VAR

    capT<-nrow(y)               # Number of observations we
                                # are going to use.

    k<-m*nsteps;                # k = mh, the maximal number

    q <- nrow(yconst)             # Number of constraints

    # Make arrays to hold the Gibbs sampler results
    yforc<-array(0,c(gibbs,nsteps,m))


    # Do the Gibbs draws.....
    for(i in 1:(burnin+gibbs))
    {
      # Step (a): Compute a draw of the conditional forecasts
      # COMPUTE INNOVATIONS: These are the structural innovations
      # Solve the constraint equation for the updated forecast errors

      # Generate the forecasts without shocks
      ytmp<-as.matrix(coef.forecast.VAR(y, intercept, ar.coefs, exog.coefs, m, p,
                                     capT, nsteps, A0=(A0)))

      # Get the impulse responses that correspond to the forecasted data
      M <- irf.VAR(varobj, nsteps, A0)$mhat

      # Construct the draw of the orthogonalized innovations that
      # satisfy the hard condition.

      # These are the q constrained innovations
      r <- (yconst - ytmp[(capT+1):(capT+nsteps),])
      r<-matrix(r[1:nsteps,1],ncol=1)

      # Build the matrix of the impulses that define the constraint
      R <- matrix(0, k, q)

      # Put the g'th column of the impulse into the constraint matrix,
      # such that R * epsilon = r
      for (g in 1:q)
      {
        if(g==1)
          { R[1:length(M[,1,1]), 1] <- M[,1,1] }
        else
          {
            R[,g] <- c(M[,1,g], R[,g-1])[1:k]
          }
      }

      # Solve the minimization problem for the mean and variance of
      # the constrained innovations.

      RRinv<-solve(crossprod(R))
      mean.epsilon <- R%*%RRinv%*%r;
      var.epsilon <- diag(1,nrow=k) - (R%*%RRinv%*%t(R));

      # Draw from the singular MVN pdf of the constrained innovations.

      epsilon.i <- matrix(rmultnorm(1, mean.epsilon, var.epsilon),
                          nrow=nsteps, ncol=m, byrow=T)

      # Add the innovations to the forecasts
      ytmp[(capT+1):(capT +nsteps),] <- ytmp[(capT+1):(capT +nsteps),]+epsilon.i%*%A0

      # Store forecasts that are past the burnin point
      if(i>burnin)
        { for(j in 1:m)
            { yforc[(i-burnin),(1:nsteps),j]<-ytmp[(capT+1):(capT+nsteps),j] }
        }


      # Step (b): Compute the mode of the posterior for the
      # conditional forecast distribution.  This is the "extended"
      # dataset that includes the i'th Gibbs sample forecast


      # Build the augmented LHS and RHS matrices
      # 1) Get the nsteps+p observations we need to build the lagged
      # endogenous variables for the augmented system.

      Y.update <- ytmp[(capT-p+1):nrow(ytmp),]

      # 2) Build the updated X -- this is hard because we need to get
      # the RHS lags correct.  We do this by padding the existing Y
      # and then building the lags.  This reuses the code for the lag
      # construction in the szbvar code

      X.update <- matrix(0, nsteps, ncoef)
      X.update[,ncoef] <- matrix(1, nsteps, 1)

      # Note that constant is put last here when the lags are made
      for(j in 1:p)
        {
          X.update[1:nsteps,(m*(j-1)+1):(m*j)] <- matrix(Y.update[(p+1-j):(nsteps+p-j),],ncol=m)
        }

      # Put on the exogenous regressors and make the constant the
      # first exog regressor after the AR coefs
      if(is.null(exog)==F)
        {
          X.update<-cbind(X.update,exog);
        }

      # Now, stack the original Y data and the augmented data.
      Y.update <- rbind(Y, ytmp[(capT+1):nrow(ytmp),])

      # Set up crossproducts and inverses we need
      X.update <- rbind(X, X.update)
      XX.update <- crossprod(X.update)     # Cross product of RHS variables
      hstar.update <- H0 + XX.update       # Prior + Cross product

      # Updated Regression estimates, Beta | Sigma
      B.update<-solve((hstar.update),(crossprod(X.update,Y.update) + H0[,1:m]))

      # Posterior mean of Sigma | Beta
      S.update <- (S0 + crossprod(Y.update)
                   + H0[1:m,1:m] - t(B.update)%*%(hstar.update)%*%(B.update))/(capT+nsteps+nu-m-1)

      # Posterior variance of B: VBh = diag(Sh.*.(inv((H0 + x'x)))
      hstarinv <- solve(hstar.update)
      vcv.Bh <- kronecker(S.update,hstarinv)

      # Draw from the conditional posterior pdfs of the parameters

      # This is only valid for just-identified models.
      df <- capT - m*p - m - 1 + nsteps
      wisharts <- rwishart(1, df, diag(m))

      # Generate the draws from the Wishart and the Beta
      # Wishart draw
      Sigmat <- (chol(S.update))
      Sigma.Draw <- t(Sigmat)%*%(df*solve(matrix(wisharts,m,m)))%*%Sigmat
      sqrtwish <- t(chol(Sigma.Draw))
      # Covariance of beta
      bcoefs.covar <- t(chol(vcv.Bh))

      # Draw of beta|Sigma ~ MVN(B.update, S.Update .*. Hstarinv)
      aplus <- matrix(B.update, ncol=1) +
          bcoefs.covar%*%matrix(rnorm(nrow(bcoefs.covar)), ncol=1)

      # Reshape and extract the coefs
      aplus <- matrix(aplus, ncol=m)
      aplus.coefs<-t(aplus[1:(m*p),]);          # extract the ar coefficients
      dim(aplus.coefs)<-c(m,m,p)                # push ar coefs into M x M x P array
      aplus.coefs<-aperm(aplus.coefs,c(2,1,3))  # reorder array

      intercept <- aplus[(m*p+1),]      # get drawn intercept....
      ar.coefs<-aplus.coefs             # AR coefs
      A0 <- sqrtwish

#      exog.coefs <- aplus[(m*p+2):nrow(aplus),]

      # Need to add something here to deal with the exogenous
      # regressors!


      # Print some intermediate results to capture progress....
      # and tell us that things are still running
      if (i%%1000==0)
        { cat("Gibbs Iteration = ", i, "     \n");
          if(i<=burnin)
            { cat("(Still a burn-in draw.)\n");
            }
        }
      # Back to the top of the Gibbs loop....
    }
  endtime<-date()
  # Print time stamp so we know how long everything took.
  cat("Start time : ", starttime, "\n");
  cat("End time   : ", endtime, "\n");
    # Returns a list object
    output <- list(forecast=yforc, orig.y=y) #llf=ts(llf),hyperp=c(mu,prior)))
    attr(output, "class") <- c("forecast.VAR")
    return(output)
}


# Forecasting function for MSBVAR models.
#
# 20110113 : Initial version
# 20110628 : Clean out remaining bugs on DF corrections for
#            observations
# 20120207 : Updated to work with revisions to MSBVAR

###### MSBVARfcast #####
# Forecast for one draw -- this is an internal function that only
# returns the forecasts averaged over the regimes.
#
# Inputs :
# y = data modeled + forecast steps
# k = number of forecast steps
# h = number of regimesq
# ar.coefs = array of c(m,m,p,h) for the MSBVAR
# intercepts = array of c(m,h) for the intercepts
# Sigma = array of c(m.m,h) for the variances
# shock = array of c(k,m,h) for the shocks
# ss = matrix of c(k,h) for the regimes

MSBVARfcast <- function(y, k, h, ar.coefs, intercepts, shocks, ss)
{
    # Get constants
    dims <- dim(ar.coefs)
    m <- dims[1]
    p <- dims[3]
    capT <- nrow(y)

    # Set up the object to hold the forecasts
    # Zero out forecast periods so we can cumulate sums over time,
    # lags, and regimes to get the correct weighted forecast.

    forcs <- rbind(y, matrix(0, k, m))

    for(j in 1:k)     # Loop over forecast steps
    {
        for(i in 1:h) # Now over the regime weights, multiplying by
                      # the regime weights and aggregating
        {
            forcs[capT + j,] <- forcs[capT+j,] +
                                ss[(capT+j),i]*(forcs[capT + j - 1,] %*%
                                                ar.coefs[,,1,i] + intercepts[,i] + shocks[j,,i])
            # Now add in the other lags -- be sure they are from the
            # right state!
            if(p>1) { for(lg in 2:p)
                      { forcs[capT + j, ] <- (forcs[capT + j, ] +
                                              ss[(capT+j-lg),i]*(forcs[capT + j - lg, ] %*% ar.coefs[,,lg,i]))
                    }
                  }
        }
    }

    return(forcs[(capT+1):(capT+k),])
}

# Betai2coefs
# Converts the draw Betai into arrays of AR coefficents for each
# regime and the intercepts
Betai2coefs <- function(Betai, m, p, h)
{
    ar.coefsi <- array(0, c(m,m,p,h))
    for(i in 1:h)
    {
        tmp <- t(Betai[1:(m*p),,i])
        dim(tmp) <- c(m,m,p)
        ar.coefsi[,,,i] <- aperm(tmp, c(2,1,3))  # lags then regimes.
#        ar.coefs[,,,i] <- tmp
}
    intercepts <- Betai[(m*p+1),,]
    return(list(ar.coefsi=ar.coefsi, intercepts=intercepts))
}


updateinit <- function(Y, init.model)
{
    n<-nrow(Y);	 	# of observations in data set
    m<-ncol(Y)	# of variables in data set
    p <- init.model$p

    # Compute the number of coefficients
    ncoef<-(m*p)+1;  # AR coefficients plus intercept in each RF equation
    ndum<-m+1;                # of dummy observations
    capT<-n-p+ndum;   	      # # of observations used in mixed estimation (Was T)
    Ts<-n-p;  		      # # of actual data observations @

    # Declare the endoegnous variables as a matrix.
    dat<-as.matrix(Y);

    # Create data matrix including dummy observations
    # X: Tx(m*p+1)
    # Y: Txm, ndum=m+1 prior dummy obs (sums of coeffs and coint).
	# mean of the first p data values = initial conditions
    if(p==1)
      { datint <- as.vector(dat[1,]) }
    else
      { datint<-as.vector(apply(dat[1:p,],2,mean)) }

    # Y and X matrices with m+1 initial dummy observations
    X<-matrix(0, nrow=capT, ncol=ncoef)
    Y<-matrix(0, nrow=capT, ncol=m)
    const<-matrix(1, nrow=capT)
    const[1:m,]<-0;	         # no constant for the first m periods
    X[,ncoef]<-const;

    # Build the dummy observations we need
    for(i in 1:m)
      {
        Y[ndum,i]<-datint[i];
        Y[i,i]<-datint[i];
        for(j in 1:p)
          { X[ndum,m*(j-1)+i]<-datint[i];
            X[i,m*(j-1)+i]<-datint[i];
          }
      }

    # Note that constant is put last here when the lags are made
    for(i in 1:p)
      { X[(ndum+1):capT,(m*(i-1)+1):(m*i)]<-matrix(dat[(p+1-i):(n-i),],ncol=m)
      }

     # Put on the exogenous regressors and make the constant the
     # first exog regressor after the AR coefs
    ## if(is.null(z)==F)
    ##   {
    ##     pad.z <- matrix(0,nrow=capT,ncol=ncol(z))
    ##     pad.z[(ndum+1):capT,] <- matrix(z[(p+1):n,], ncol=ncol(z))
    ##     X<-cbind(X,pad.z);
    ##   }

    # Get the corresponding values of Y
    Y[(ndum+1):capT,]<-matrix(dat[(p+1):n,],ncol=m);

    # Weight dummy observations
    X[1:m,] <- init.model$prior[6]*X[1:m,];
    Y[1:m,] <- init.model$prior[6]*Y[1:m,];
    X[ndum,] <- init.model$prior[7]*X[ndum,];
    Y[ndum,] <- init.model$prior[7]*Y[ndum,];

    init.model$X <- X
    init.model$Y <- Y
    init.model$hstar <- (init.model$H0 + crossprod(X))

    return(init.model)
}


# This is an amalgam of uc.forecast and gibbs.msbvar.  In this
# version, the same steps as gibbs.msbvar are used, but with the
# updated data from the forecast step.

# x = posterior mode object -- can either be from gibbs.msbvar or
#     msbvar.  Need to handle identification steps or permute as
#     appropriate.  Should warn against permuting!
# k = number of forecast steps.
#
# -----------------------------------------------------------------#
################## Steps in MSBVAR forecasting #####################
# -----------------------------------------------------------------#
# A) Augment the dataset from the mode -- forecasting
# B) Update the input data
# C) For given parameters,
#    1) draw statespace
#    2) update / draw Q
#    3) update regression -- need to remake the input matrices as part
#       of this.
#    4) sample variances
#    5) sample regression
#    6) update errors
#    7) Impose identification on the regimes (if necessary)
# top of loop
# -----------------------------------------------------------------#
#
# Initial version assumes user is inputing an identified object from
# gibbs.msbvar().  Need to add a warning when they using an msbvar()
#     object and reconcile this.
#
forecast.MSBVAR <- function(x, k, N1=1000, N2=1000)
{
    # Get some constants / objects
    init.model <- x$init.model  # initial model with the prior
                                # matrices.

    # Get other constants from the input object
    h <- x$h
    m <- x$m
    p <- x$p

    # Set the sampler method for Q based on inputs
    if(attr(x, "Qsampler")=="Gibbs") Qsampler <- Q.drawGibbs
    if(attr(x, "Qsampler")=="MH") Qsampler <- Q.drawMH

    # Input data
    y <- init.model$y

    # Get sample size of original dataset and the effective sample
    TT <- capT <- nrow(init.model$y)

    # Settle how we are going to handle input sample size versus the
    # other sample sizes for the inputs / outputs.  Do this the same
    # way as the earlier unconditional forecasting code.  In this
    # code, we treated capT as the full original sample size minus the
    # lags to make it match the estimator. We then
    # manipulate from that.  We do this by always forecasting off of
    # the original data from capT+1 to capT+k,  This is the principal
    # in uc.forecast() and forecast.VAR and coef.forecast.VAR() [in
    # hidden.R].  The difference here is that tbe value of capT has to
    # be adjusted compared to what is in the other code.
    #
    # Principal in this is just working with data augmentation off of
    # the original data setup.

    Tp <- TT-p # effective sample

    alpha.prior <- x$alpha.prior

    # Get the modes for the posterior
    Q <- matrix(apply(x$Q.sample, 2, mean), h, h)
    Betai <- array(apply(x$Beta.sample, 2, mean), c(m*p+1, m, h))
    tmp <- Betai2coefs(Betai, m, p, h)
    ar.coefsi <- tmp$ar.coefsi
    intercepts <- tmp$intercepts

    Sigmai <- array(apply(matrix(apply(x$Sigma.sample, 2, mean),
                                 m*(m+1)/2, h),
                          2, xpnd), c(m,m,h))
    intercepts <- Betai[(m*p+1),,]

    # Extend the dataset / input data
    # Residuals and regression
    hreg <- hregime.reg2(h, m, p, mean.SS(x), init.model)
    e <- array(0, c(Tp+k, m, h))

    # Now take random draws from the error process for each regime to
    # pad out initial values for the forecast periods.
    Sigmachol <- aperm(array(apply(Sigmai, 3, chol), c(m,m,h)),
                       c(2,1,3))

    for(i in 1:h)
    {
        etmp <-  t(Sigmachol[,,i]%*%matrix(rnorm(k*m), m, k))
        e[,,i] <- rbind(hreg$e[,,i], etmp)
    }

    # State-space initialization -- flat convolution of the last
    # state!

    tmpSS <- mean.SS(x)
    sstmp <- rbind(mean.SS(x), matrix(rep(NA, k*h), k, h))

    for(i in 1:k)
    {
        tmp <- sstmp[Tp+i-1,]%*%Q
        sstmp[Tp+i,] <- tmp/sum(tmp)

    }

    transtmp <- count.transitions(sstmp)
    ss <- list(SS=sstmp, transitions=transtmp)

    # Make list objects for Beta and Sigma
    Betai <- list(Betai=Betai)
    Sigmai <- list(Sigmai=Sigmai)

    # Storage
    forecasts <- array(0, c(N2, k, m))
    states <- vector("list", N2)

    # Burnin loop + Final loop
    for (j in 1:(N1+N2))
    {

        # Update the residuals
        shocks <- array(rnorm(m*k*h), c(k, m, h))

        Sigmai.chol <- aperm(array(apply(Sigmai$Sigmai, 3, chol), c(m,m,h)),
                             c(2,1,3))

        for(i in 1:h)
         {
#             shocks[,,i] <- Sigmai.chol[,,i]%*%shocks[,,i]
             shocks[,,i] <- t(tcrossprod(Sigmai.chol[,,i], shocks[,,i]))
             e[(Tp+1):(Tp+k),,i] <- (shocks[,,i])
         }


        # Forecast

        fcast <- MSBVARfcast(y[(p+1):nrow(y),], k, h, ar.coefsi,
                             intercepts, shocks, ss$SS)

#        print(round(cbind(shocks[,,1], shocks[,,2], fcast), 2))

        # Update / Draw statespace

        oldtran <- matrix(0,h,h)

        while(sum(diag(oldtran)==0)>0)
        {
            ss <- SS.ffbs(e, (Tp+k), m, p, h, Sigmai$Sigmai, Q)
            oldtran <- ss$transitions
        }

        #print(ss$transitions)
        ## plot(ts(ss$SS), plot.type="s", col=1:h,
        ##      main=paste("Iteration :", j))

        #if(sum(diag(ss$transitions)==0)>0) stop("Oops.")

        # Draw Q
        Q <- Qsampler(Q, ss$transitions, prior=alpha.prior, h)

        # Update regression
        # Update the model object input with the new data -- this is
        # not very efficient right now, but it works

        # Need to re-initialize the init.model object with the latest
        # forecast data so it can be an input to the next sample calls
        #
        # This should be really done with some BVAR setup function
        # that is then used in all of the later models in the pkg --
        # something to do later for speed.

        init.model <- updateinit(Y=rbind(y,fcast), init.model)

        # Update regression steps
        hreg <- hregime.reg2(h, m, p, ss$SS, init.model)

        # Draw variances

#        cat("Iteration :", j, "\n")
#        print(hreg$Sigmak)

        Sigmai <- Sigma.draw(m, h, ss, hreg, Sigmai$Sigmai)

#        print(Sigmai)

        # Draw VAR regression coefficients
        Betai <- Beta.draw(m, p, h, Sigmai$Sigmai, hreg, init.model,
                           Betai$Betai)

        # Split out AR coefs and intercepts
        tmp <- Betai2coefs(Betai$Betai, m, p, h)
        ar.coefsi <- tmp$ar.coefsi
        intercepts <- tmp$intercepts

        # Update error estimate
        e <- residual.update(m, h, init.model, Betai$Betai, e)

        # Permute regimes if necessary / order regimes

        # Save results if past burnin=N1.
        # Store forecasts and states
        if(j > N1)
        {
            for(i in 1:m)
            {
                forecasts[(j-N1),,i] <- t(fcast[,i])
                states[[(j-N1)]] <-
                    as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
            }
        }

        # Print out iteration information
        if(j%%1000==0) cat("Iteration : ", j, "\n")
    }

    # Define the output object and its attributes
    class(forecasts) <- c("forecast.MSBVAR")
    class(states) <- c("SS")

    output <- list(forecasts=forecasts, ss.sample=states, k=k, h=h)

    # Add classing and attributes here.  These will be use later for
    # the plotting and other summary functions.
    attr(output, "eqnames") <- attr(x, "eqnames")
    tmp <- tsp(x$init.model$y)
    attr(output, "start") <- tmp[1]
    attr(output, "end") <- tmp[2]
    attr(output, "freq") <- tmp[3]

    return(output)
}

