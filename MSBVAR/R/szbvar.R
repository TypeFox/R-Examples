
# szbvar estimator

"szbvar" <-
function(Y, p, z=NULL, lambda0, lambda1, lambda3,
                 lambda4, lambda5, mu5, mu6, nu=ncol(Y)+1, qm=4, prior=0,
                 posterior.fit=FALSE)
{

    sanity.check.bvar(list(Y=Y, p=p, z=z, lambda0=lambda0,
                           lambda1=lambda1, lambda3=lambda3,
                           lambda4=lambda4, lambda5=lambda5,
                           mu5=mu5, mu6=mu6, qm=qm, prior=prior,
                           posterior.fit=posterior.fit))

    n<-nrow(Y);	 	# # of observations in data set
    m<-ncol(Y);	 	# # of variables in data set

    # Test for the exogenous variables
    if (is.null(z))
      { num.exog <- 0 }
    else
      { num.exog <- ncol(z)
        z <- as.matrix(z)
      }

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
      {   # wrap first p values in an as.matrix() in case someone
          # passes in a univariate series.
          datint<-as.vector(apply(as.matrix(dat[1:p,]),2,mean))
      }

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
    if(is.null(z)==F)
      {
        pad.z <- matrix(0,nrow=capT,ncol=ncol(z))
        pad.z[(ndum+1):capT,] <- matrix(z[(p+1):n,], ncol=ncol(z))
        X<-cbind(X,pad.z);
      }

    # Get the corresponding values of Y
    Y[(ndum+1):capT,]<-matrix(dat[(p+1):n,],ncol=m);

    # Weight dummy observations
    X[1:m,]<-mu5*X[1:m,];
    Y[1:m,]<-mu5*Y[1:m,];
    X[ndum,]<-mu6*X[ndum,];
    Y[ndum,]<-mu6*Y[ndum,];

    # END OF DATA SET UP FOR THE MIXED ESTIMATION

    # NOW CREATE THE PRIORS

    # Create monthly lag decay to match 1/k decay in quarterly
    # data where k = quarters
    ld<-seq(1:p)^-lambda3;			# regular lag decay (note ^-lambda3)
    if(qm==12)
      { j<-ceiling(p/3)^-lambda3;   	# last quarter (rounded up) eg. l2=13=>xx2=5
        b<-0
        if(p > 1)
          { b<-(log(1)-log(j))/(1-p) }
        a<-exp(-b);
        ld<-a*exp(b*seq(1:p));	# Tao Zha's lag decay to match 13th lag
      }

    # Scale factors from OLS
    s2<-matrix(0,nrow=m,ncol=1)
    for(i in 1:m)
      {
        s2[i,1] <- ar.ols(Y[,i], aic=FALSE, order.max=p,
                          intercept=TRUE, demean=FALSE)$var.pred
      }
    #  Prior scale matrix for Sigma:  Sigma ~ IW(H0,v) @
    S0 <- diag(m)
    diag(S0)<- s2/(lambda0^2);	 # see SZ p.961 (24)

    # prior for intercept
    prior.intercept <- 0;
    if( lambda4 > 0)
      { prior.intercept <- 1/(lambda0*lambda4)^2 }  # inverse prior variance for constant @

    # Now scale the priors for the exogneous variables
    prior.exog <- 0
    if( lambda5 > 0 )
      { prior.exog <- 1/(lambda0*lambda5)^2 }


    # Prior cov of B|Sigma is Sigma.*.inv(H0)  for now, this uses the
    # prior on the intercept for the exogenous variables

    if(num.exog==0)
      {  H0<-diag(c(kronecker((1/(ld*lambda0*lambda1))^2,s2),
                    prior.intercept), nrow=ncoef, ncol=ncoef)
       }
    else
      {  H0<-diag(c(kronecker((1/(ld*lambda0*lambda1))^2,s2),
                    prior.intercept, rep(prior.exog, num.exog)),
                  nrow=(ncoef+num.exog), ncol=(ncoef+num.exog))
       }

    # Now, do the special cases of normal-flat and flat-flat
    if(prior == 1)
      { S0 <-0*S0 } # using normal-flat prior as special case

    if(prior == 2)         # using the flat-flat prior as a special case
      { H0<-0*H0
        S0<-0*S0
        nu <- 0
      }

    # Set up some matrices for later
    XX <- crossprod(X)    # Cross product of RHS variables
    hstar1 <- H0 + XX     # Prior + Cross product

    # Posterior mean of B: Bh  = inv(x'x + H0)*(x'y + H0[.,1:m]);
    # This is different from the original code because the intercept
    # is now the last coefficient in the matrix, rather than the
    # first, followed by the exogenous regressors.
    Bh<-solve((hstar1),(crossprod(X,Y) + H0[,1:m]))

    # Posterior mean of Sigma
    # This is different from the original code because the intercept
    # is now the last coefficient in the matrix, rather than the first.
    St <- (S0 + crossprod(Y) + H0[1:m,1:m] - t(Bh)%*%(hstar1)%*%(Bh))
    Sh<- St/(Ts+nu-m-1)

    # Posterior variance of B: VBh = diag(Sh.*.(inv((H0 + x'x)))
    hstarinv <- solve(hstar1)
    vcv.Bh <- kronecker(Sh,hstarinv)

    # Residuals and MLE estimators of the variance
    u<-(Y - X%*%(Bh));
    Sh1<-crossprod(u)/capT;  # u'u/n form of estimator

    # Format the output variables
    # Split the coefficient matrix (will make IRFs and forecasting easier)
    intercept<-Bh[ncoef,];                   # extract the intercept
    ar.coefs<-t(Bh[1:(ncoef-1),]);           # extract the
    # ar coefficients

    dim(ar.coefs)<-c(m,m,p)              # push ar coefs into M x M x P array

    ar.coefs<-aperm(ar.coefs,c(2,1,3))   # reorder array so columns are for eqn

    if(is.null(z))
      {exog.coefs <- NA}
    else
      {
        exog.coefs <- Bh[(ncoef+1):nrow(Bh),]
      }

    marg.llf <- NA
    marg.post <- NA
    coef.post <- NA

    pfit <- list(capT=capT, m=m, ncoef=ncoef, num.exog=num.exog,
                          nu=nu, H0=H0, S0=S0, Y=Y, X=X,
                          hstar1=hstar1, Sh=Sh, u=u, Bh=Bh, Sh1=Sh1)

    output <- list(intercept = intercept,
                   ar.coefs=ar.coefs,
                   exog.coefs=exog.coefs,
                   Bhat = Bh,
                   vcv=Sh1,
                   vcv.Bh=vcv.Bh,
                   mean.S=Sh,
                   St=St,
                   hstar=(H0 + XX),
                   hstarinv=hstarinv,
                   H0=H0,
                   S0=S0,
                   residuals = u,
                   X=X,
                   Y=Y,
                   y=dat,
                   z=z,
                   p=p,
                   num.exog=num.exog,
                   qm=qm,
                   prior.type=prior,
                   prior=c(lambda0,lambda1,lambda3,lambda4,lambda5,mu5,mu6,nu),
                   pfit=pfit,
                   marg.llf=marg.llf,
                   marg.post=marg.post,
                   coef.post=coef.post)
    class(output) <- c("BVAR")
    attr(output, "eqnames") <- colnames(dat) # Get variable names for
                                             # attr

    # Compute the posterior fit measures.
    if(posterior.fit==T)
    {
        tmp <- posterior.fit.BVAR(output)
        output$marg.llf <- tmp$data.marg.llf
        output$marg.post <- tmp$data.marg.post
        output$coef.post <- tmp$coef.post
    }

    # Here are the returns
    return(output)
  }


# Summary function for BVAR models
"summary.BVAR" <- function(object, ...)
{
    cat("------------------------------------------\n")
    cat("Sims-Zha Prior reduced form Bayesian VAR\n")
    cat("------------------------------------------\n")
    if(object$prior.type==0) prior.text <- "Normal-inverse Wishart"
    if(object$prior.type==1) prior.text <- "Normal-flat"
    if(object$prior.type==2) prior.text <- "Flat-flat"

    cat("Prior form : ", prior.text, "\n")
    cat("Prior hyperparameters : \n")
    cat("lambda0 =", object$prior[1], "\n")
    cat("lambda1 =", object$prior[2], "\n")
    cat("lambda3 =", object$prior[3], "\n")
    cat("lambda4 =", object$prior[4], "\n")
    cat("lambda5 =", object$prior[5], "\n")
    cat("mu5     =", object$prior[6], "\n")
    cat("mu6     =", object$prior[7], "\n")
    cat("nu      =", object$prior[8], "\n")

    cat("------------------------------------------\n")
    cat("Number of observations : ", nrow(object$Y), "\n")
    cat("Degrees of freedom per equation : ", nrow(object$Y)-nrow(object$Bhat), "\n")
    cat("------------------------------------------\n")

    cat("Posterior Regression Coefficients :\n")
    cat("------------------------------------------\n")
    cat("Autoregressive matrices: \n")
    for (i in 1:dim(object$ar.coefs)[3])
    {
        cat("B(", i, ")\n", sep="")
        prmatrix(round(object$ar.coefs[,,i], 6))
        cat("\n")
    }
    cat("------------------------------------------\n")
    cat("Constants\n")
    cat(round(object$intercept,6), "\n")
    cat("------------------------------------------\n")

    if(is.na(object$exog.coefs[1])==FALSE)
    {
        cat("------------------------------------------\n")
        cat("Exogenous variable posterior coefficients\n")
        prmatrix(object$exog.coefs)
        cat("\n")
        cat("------------------------------------------\n")
    }

    cat("------------------------------------------\n")
    cat("Posterior error covariance\n")
    prmatrix(object$mean.S)
    cat("\n")
    cat("------------------------------------------\n")

}
