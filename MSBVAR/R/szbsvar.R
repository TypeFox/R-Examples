"szbsvar" <-
function(Y, p, z=NULL, lambda0, lambda1, lambda3, lambda4, lambda5,
         mu5, mu6, ident, qm=4){

    sanity.check.bsvar(list(Y=Y, p=p, z=z, lambda0=lambda0,
                           lambda1=lambda1, lambda3=lambda3,
                           lambda4=lambda4, lambda5=lambda5,
                           mu5=mu5, mu6=mu6, qm=qm, ident=ident))

    # Set up some constants we need
    m <- ncol(Y)                                  # number of endog variables
    nexog <- ifelse(is.null(z)==TRUE, 0, ncol(z)) # number of exog variables
    ncoef <- m*p + nexog + 1                      # plus 1 for the constant
    n<-nrow(Y);	 	                    # observations

    if(dim(ident)[1]!=m)
    {
        stop("Identification matrix 'ident' and dimension of 'Y' are nonconformable.")
    }


      # Compute the number of coefficients
    endog.ncoef<-(m*p)+1;  # AR coefficients plus intercept in each RF equation
    ndum<-m+1;             # of dummy observations
    capT<-n-p+ndum;        # degrees of freedom for the mixed estimation
    Ts<-n-p;  	     # # of actual data observations

      # Do some error checking of the identification scheme.
      # NEED to put this into the model.....

      # Define the linear restrictions from the identification matrix
      # for the free parameters

      # set up the linear restrictions for the parameters based on the A0
      # identification in ident
    Q <- array(0, c(m,m,m))
    for(i in 1:m)
    { Q[,,i] <- diag(ident[,i])
  }

      # Find the orthonormal bases for each of the Q matrices.  Need
      # thse because they define the null space for the squeezing the
      # parameters. This is just the null space for each matrix in Q.

    Ui <- sapply(1:m, function(i){null.space(Q[,,i])}, simplify=F)

      # Set up the prior

      # Prior mean of the regression parameters, Pi, for the ith equation
    Pi <- matrix(0, ncoef, m)
    diag(Pi) <- 1

      # Set up the prior for each equation using the resids from a
      # univariate AR model

      # Scale factors from univariate OLS
    s2i<-matrix(0,nrow=m,ncol=1)
    for(i in 1:m)
    {
        s2i[i,1] <- ar.ols(Y[,i], aic=FALSE,order.max=p,
                           intercept=TRUE,demean=FALSE)$var.pred
    }

      # Prior scale for A0 -- Si in the Waggoner and Zha notation
    S0 <- diag(m)
    diag(S0)<- 1/s2i

      # Prior for A+ or F coefficients -- same for all equations for now

      # Lag decays
      # create monthly lag decay to match 1/k decay in quarterly
      #  data where k = quarters
    if(qm==12)
    { j<-ceiling(p/3)^-lambda3;   	# last quarter (rounded up) eg. l2=13=>xx2=5
      b<-0
      if(p > 1)
      { b<-(log(1)-log(j))/(1-p) }
      a<-exp(-b);
  }

      # Find the lag decays for the variances over the p lags
    Aplus.prior.cov <- matrix(0, ncoef, 1)

      for(i in 1:p)
        { if(qm==12)
            { ld <- a*exp(b*i*lambda3) }
          for(j in 1:m)
            { if(qm==12)
                { Aplus.prior.cov[((i-1)*m+j),1] <- ld^2/s2i[j,1]
                } else {
                  Aplus.prior.cov[((i-1)*m+j),1] <- (1/i^lambda3)^2/s2i[j,1]
                }
            }
        }

      # Find the prior for A0 conditional prior variances
      A0.prior <- lambda0^2/s2i           # sg0bida
      Aplus.prior <- lambda0^2*lambda1^2*Aplus.prior.cov #sgpbida
      Aplus.prior[(m*p+1),1] <- (lambda0*lambda4)^2  # prior for intercept

      if(nexog>0)
        { Aplus.prior[(m*p+2):ncoef,1] <- lambda0^2*lambda5^2  # prior for eexog
        }

      Aplus.prior1 <- Aplus.prior[(m+1):ncoef,1] # sgppbd

      # Now compute the H matrices
      Hptd <- diag(as.vector(Aplus.prior))
      Hptdi <- diag(as.vector(1/Aplus.prior))

      # Now find the final covariance matrices for each of the i=1,..,m
      # equations.
      H0multi <- array(0, c(m, m, m))
      H0invmulti <- H0multi
      Hpmulti <- array(0, c(ncoef,ncoef,m))
      Hpmultiinv <- Hpmulti

      # This can be modified if we want to implement an asymmetric prior
      # across the columns (equations).
      H0td <- matrix(0, m, m)
      H0tdi <- H0td
      for (i in 1:m)
        { # A0 parts
          A0i <- A0.prior
          A0i.inv <- 1/A0i
          diag(H0td) <- A0i
          diag(H0tdi) <- 1/A0i
          H0multi[,,i] <- H0td
          H0invmulti[,,i] <- H0tdi

          # A+ parts
          Hpmulti[,,i] <- Hptd
          Hpmultiinv[,,i] <- Hptdi
        }

      # Now combine the prior with the linear restrictions --- maps from
      # a(i) and f(i) vectors into the b(i) and g(i) vectors of the
      # restricted model.  This is where we map from q(a_i, f_i) to
      # q(a_i, f_i | Q a_i = 0 and R f_i = 0)

      # This is where we make the Ptilde, H0tilde and Hptilde matrices for
      # the restricted system.

      Hpinv.tilde <- Hpmultiinv

      # Use sapply, since we do not know the size of the null spaces, so the
      # result needs to be dynamically sized

      Pi.tilde <- sapply(1:m, function(i){Pi%*%Ui[[i]]}, simplify=F)
      H0inv.tilde <- sapply(1:m, function(i)
                            {t(Ui[[i]])%*%H0invmulti[,,i]%*%Ui[[i]]}, simplify=F)

      # Set up the data
      # 1) Set up matrices of data and the moment matrices for the data

      # Test for the exogenous variables and check rank
      if (is.null(z))
        { num.exog <- 0 } else {
            num.exog <- ncol(z)
            z <- as.matrix(z)
            if(det(crossprod(cbind(rep(1,nrow(z)),z)))<=0)
            {
                stop("Matrix of exogenous variables, z has deficient rank.")
            }
        }

      # Create data matrix including dummy observations
      # X: Tx(m*p+1)
      # Y: Txm, ndum=m+1 prior dummy obs (sums of coeffs and coint).
      # mean of the first p data values = initial conditions
      if(p==1)
        { datint <- as.vector(Y[1,])
        } else {
          datint<-as.vector(apply(Y[1:p,],2,mean)) }

      # Y and X matrices with m+1 initial dummy observations
      X1<-matrix(0, nrow=capT, ncol=endog.ncoef)
      Y1<-matrix(0, nrow=capT, ncol=m)
      const<-matrix(1, nrow=capT)
      const[1:m,]<-0;	         # no constant for the first m periods
      X1[,endog.ncoef]<-const;

      # Build the dummy observations we need
      for(i in 1:m)
        {
          Y1[ndum,i]<-datint[i];
          Y1[i,i]<-datint[i];
          for(j in 1:p)
            { X1[ndum,m*(j-1)+i]<-datint[i];
              X1[i,m*(j-1)+i]<-datint[i];
            }
        }

      # Make the lags.  Note that constant is put last here when the lags are made
      for(i in 1:p)
        { X1[(ndum+1):capT,(m*(i-1)+1):(m*i)]<-matrix(Y[(p+1-i):(n-i),],ncol=m)
        }

      # Put on the exogenous regressors and make the constant the
      # first exog regressor after the AR coefs
      if(is.null(z)==F)
        {
          pad.z <- matrix(0,nrow=capT,ncol=ncol(z))
          pad.z[(ndum+1):capT,] <- matrix(z[(p+1):n,], ncol=ncol(z))
          X1<-cbind(X1,pad.z);
        }

      # 2) Dummy observations

      # Get the corresponding values of Y
      Y1[(ndum+1):capT,]<-matrix(Y[(p+1):n,],ncol=m);

      # Weight dummy observations
      X1[1:m,]<-mu5*X1[1:m,];
      Y1[1:m,]<-mu5*Y1[1:m,];
      X1[ndum,]<-mu6*X1[ndum,];
      Y1[ndum,]<-mu6*Y1[ndum,];

      # 3) Get moment matrices
      XX <- crossprod(X1)
      XY <- crossprod(X1, Y1)
      YY <- crossprod(Y1)

      # Compute the posterior moments -- combine the data and the
      # prior.  These are the posterior moments needed for the results
      # in Waggoner and Zha 2003, JEDC, eqns. 12 and 13 (the H_iu,
      # P_i, S_i below these equations

      # H_i^-1 matrices for equations 1...m
      Hpinv.posterior <- sapply(1:m, function(i) { XX + Hpinv.tilde[,,i] },
                                simplify=F)

      # Next two build the P_i, the "squeezed" matrices of the SVAR parameters
      P1.posterior <- sapply(1:m,
                             function(i) { XY%*%Ui[[i]] +
                                             Hpinv.tilde[,,i]%*%Pi.tilde[[i]] }, simplify=F)

      P.posterior <- sapply(1:m, function(i)
                            {solve(Hpinv.posterior[[i]])%*%P1.posterior[[i]] },
                            simplify=F)

      # S_i matrices are next -- these are the SVAR covariances for
      # the columns of A0, a_i.  Note that we do NOT include the
      # scaling by 1/T because this is not necessary for the posterior
      # draws.

      H0inv.posterior <- sapply(1:m, function(i)
                                { (t(Ui[[i]])%*%YY%*%Ui[[i]] + H0inv.tilde[[i]] +
                                   t(Pi.tilde[[i]])%*%Hpinv.tilde[,,i]%*%Pi.tilde[[i]] -
                                   t(P1.posterior[[i]])%*%P.posterior[[i]])
                          }, simplify=F)

      # Optimize the likelihood to solve for A0 parameters (the b's in
      # the "squeezed" model.  Here we are numerically computing the
      # peak of the PDF -- makes for more efficient draws from the
      # posterior later.

      # Find # free parameters in each equation and a vector of the
      # cusum for indexing them

      n0 <- sapply(1:m, function(i) {ncol(as.matrix(Ui[[i]]))})
      n0cum <- c(0,cumsum(n0))
      # Generate some random starting values.
      b <- (1/max(s2i))*(rnorm(sum(n0)))

      # Optimize the log posterior wrt A0.

      # Start with some Nelder-Mead simplex steps to get things
      # started in the right direction.
      cat("Estimating starting values for the numerical optimization\nof the log posterior of A(0)\n")

      max.obj <- optim(b, A0.llf, method=c("Nelder-Mead"),
                       control=list(maxit=6000, fnscale=capT, trace=0),
                       Ui=Ui, df=capT, H0inv.posterior=H0inv.posterior)

      # Do the final optimization with BFGS.
      cat("Estimating the final values for the numerical optimization\nof the log posterior of A(0)\n")

      max.obj <- optim(max.obj$par, A0.llf, method=c("BFGS"), hessian=F,
                       control=list(maxit=5000, fnscale=capT, trace=1),
                       Ui=Ui, df=capT,
                       H0inv.posterior=H0inv.posterior)

      # Check for convergence
      if(max.obj$convergence!=0)
        {
          stop("Estiamtes of A(0) did not converge.  You should restart the function with a new seed.")
        }


      # Build back the A0 from the b's estimated at the peak of the
      # log posterior pdf.

      # Estimate of A0 at the peak of the log-posterior.
      A0.mode <- b2a(max.obj$par, Ui)

      # Estimates of the SVAR posterior coefs for the lagged and
      # exogenous variables -- the F matrix defined in equations 13
      # and on page 351.

      # Build back the A+ and F coefficients from the sub-space of the SVAR
      # back to the unrestricted parameter space

      F.posterior <- matrix(0, ncoef, m)

      for (i in 1:m)
        { bj <- max.obj$par[(n0cum[i]+1):(n0cum[(i+1)])]
          gj <- P.posterior[[i]]%*%bj
          F.posterior[,i] <- gj
        }

      # Now map it all back to reduced form coefficients
      B.posterior <- F.posterior%*%solve(A0.mode)

      # Pluck out the arrays of the AR coefficients

      AR.coefs.posterior <- t(B.posterior[1:(m*p),])
      dim(AR.coefs.posterior) <- c(m,m,p)
      AR.coefs.posterior <- aperm(AR.coefs.posterior, c(2,1,3))

      # compute the structural innovations
      structural.innovations <- Y1%*%A0.mode - X1%*%F.posterior

      # reduced form exogenous coefficients
      if(nexog==0){
          exog.coefs <- NA
      } else {
          exog.coefs <- B.posterior[((m*p)+2):nrow(B.posterior),]
      }

    # gc() call for some cleanup
    gc(); gc();

      # Now build an output list / object for the B-SVAR model
    output <- list(XX=XX,                               # data matrix moments with dummy obs
                   XY=XY,
                   YY=YY,
                   y=Y,
                   Y=Y1,
                   X=X1,
                   structural.innovations=structural.innovations,
                   Ui=Ui,                                         # restriction transformation
                   Hpinv.tilde=Hpinv.tilde,    # Prior moments
                   H0inv.tilde=H0inv.tilde,
                   Pi.tilde=Pi.tilde,
                   Hpinv.posterior=Hpinv.posterior,
                   P.posterior=P.posterior,
                   H0inv.posterior=H0inv.posterior,
                   A0.mode=A0.mode,
                   F.posterior=F.posterior,
                   B.posterior=B.posterior,
                   ar.coefs=AR.coefs.posterior,
                   intercept=B.posterior[(m*p+1),],
                   exog.coefs=exog.coefs,
                   prior=c(lambda0,lambda1,lambda3,lambda4,lambda5,mu5,mu6),
                   df=capT,
                   n0=n0,
                   ident=ident,
                   b=max.obj$par
                   )
    class(output) <- c("BSVAR")
    attr(output, "eqnames") <- colnames(Y) # Get variable names for
                                           # attr

    return(output)
}


# Summary function for BSVAR models
"summary.BSVAR" <- function(object, ...)
{
    # Get p
    p <- dim(object$ar.coefs)[3]

    cat("------------------------------------------\n")
    cat("A0 restriction matrix\n")
    cat("------------------------------------------\n")
    prmatrix(object$ident)
    cat("\n")

    cat("------------------------------------------\n")
    cat("Sims-Zha Prior Bayesian Structural VAR\n")
    cat("------------------------------------------\n")
##     if(object$prior.type==0) prior.text <- "Normal-inverse Wishart"
##     if(object$prior.type==1) prior.text <- "Normal-flat"
##     if(object$prior.type==2) prior.text <- "Flat-flat"

    cat("Prior form : Sims-Zha\n")
    cat("Prior hyperparameters : \n")
    cat("lambda0 =", object$prior[1], "\n")
    cat("lambda1 =", object$prior[2], "\n")
    cat("lambda3 =", object$prior[3], "\n")
    cat("lambda4 =", object$prior[4], "\n")
    cat("lambda5 =", object$prior[5], "\n")
    cat("mu5     =", object$prior[6], "\n")
    cat("mu6     =", object$prior[7], "\n")
    cat("nu      =", dim(object$ar.coefs)[1]+1, "\n")

    cat("------------------------------------------\n")
    cat("Number of observations : ", nrow(object$Y), "\n")
    cat("Degrees of freedom per equation : ", nrow(object$Y)-nrow(object$Bhat), "\n")
    cat("------------------------------------------\n")

    cat("Posterior Regression Coefficients :\n")
    cat("------------------------------------------\n")
    cat("Reduced Form Autoregressive matrices: \n")
    for (i in 1:dim(object$ar.coefs)[3])
    {
        cat("B(", i, ")\n", sep="")
        prmatrix(round(object$ar.coefs[,,i], 6))
        cat("\n")
    }
    cat("------------------------------------------\n")
    cat("Reduced Form Constants\n")
    cat(round(object$intercept,6), "\n")
    cat("------------------------------------------\n")

    if(nrow(object$B.posterior)>m*p + 1)
    {
        cat("------------------------------------------\n")
        cat("Reduced Form Exogenous coefficients\n")
        prmatrix(object$B.posterior[(m*p+2):nrow(object$B.posterior),])
        cat("\n")
        cat("------------------------------------------\n")
    }

    # Now print the structural coefficients in the same way as the
    # RFs.
    cat("Structural Autoregressive matrices: \n")
    m <- dim(object$ar.coefs)[1]
    p <- dim(object$ar.coefs)[3]
    struct.ar <- object$F.posterior[1:(m*p),]
    dim(struct.ar) <- c(m, m, p)
    struct.ar <- aperm(struct.ar, c(2,1,3))
    for (i in 1:p)
    {
        cat("A(", i, ")\n", sep="")
        prmatrix(round(struct.ar[,,i], 6))
        cat("\n")
    }
    cat("------------------------------------------\n")
    cat("Structural Constants\n")
    cat(round(object$F.posterior[(m*p +1),],6), "\n")
    cat("------------------------------------------\n")

    if(nrow(object$B.posterior)>m*p + 1)
    {
        cat("------------------------------------------\n")
        cat("Structural Exogenous coefficients\n")
        prmatrix(object$F.posterior[(m*p+2):nrow(object$F.posterior),])
        cat("\n")
        cat("------------------------------------------\n")
    }


    cat("------------------------------------------\n")
    cat("Posterior mode of the A0 matrix\n")
    prmatrix(round(object$A0.mode,6))
    cat("\n")
    cat("------------------------------------------\n")

}
