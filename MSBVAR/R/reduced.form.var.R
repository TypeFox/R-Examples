"reduced.form.var" <- function(Y, p, z=NULL)
  {
      sanity.check.var(list(Y=Y, p=p, z=z))

      y <- Y
      dat<-as.matrix(Y);
      n<-nrow(dat);	 	# of observations in data set
      m<-ncol(dat);	 	# of variables in data set
      ncoef<-(m*p)+1;	        # coefficients in each RF equation
      ndum<-m+1;
      capT<-n-p+ndum;
      Ts<-n-p;  		# of actual data observations

      # Y and X matrices
      X<-matrix(0,nrow=Ts,ncol=ncoef)
      Y<-matrix(0,nrow=Ts,ncol=m)
      const<-matrix(1,nrow=Ts)
      X[,ncoef]<-const;

      # Note that constant is put last here
      for(i in 1:p)
      { X[(1:Ts),(m*(i-1)+1):(m*i)]<-matrix(dat[(p+1-i):(n-i),],ncol=m) }

      X<-cbind(X[,1:ncoef-1],X[,ncoef])

      Y[1:Ts,]<-matrix(dat[(p+1):n,],ncol=m);

      if (is.null(z)==FALSE)
        { X <- cbind(X, z[(p+1):n,]) }
      # Set up some matrices for later
      XX <- crossprod(X)
      Bh<-solve(XX, crossprod(X,Y), tol = 1e-16)
      u<-(Y - X%*%(Bh));
      Sh<-crossprod(u)/Ts;  # u'u/n form of estimator
      Sh1<-crossprod(u)/capT
      # Format the output variables
      # Split the coefficient matrix (will make IRFs and forecasting easier)
      intercept<-Bh[(m*p+1),];             # extract the intercept
      ar.coefs<-t(Bh[1:(m*p),]);           # extract the ar coefficients
      dim(ar.coefs)<-c(m,m,p)              # push ar coefs into M x M x P array
      ar.coefs<-aperm(ar.coefs,c(2,1,3))   # reorder array so columns are for eqn

      if(is.null(z)==FALSE)
        {
          exog.coefs <- Bh[(m*p+2):nrow(Bh),]
          num.exog <- ncol(z)
          z <- as.matrix(z)
        }
      else{
          exog.coefs <- NA
          num.exog <- 0
      }

      pfit <- list(capT=capT, ncoef=ncoef, num.exog=num.exog,
                   Sh1=Sh1)



      output <- list(intercept = intercept,
                     ar.coefs=ar.coefs,
                     Bhat = Bh,
                     vcv=Sh,
                     exog.coefs=exog.coefs,
                     residuals=u,
                     mean.S=Sh,
                     hstar=XX,
                     X=X,
                     Y=Y,  # the VAR Y
                     y=y, # original y
                     pfit=pfit)
      class(output) <- c("VAR")
      attr(output, "eqnames") <- colnames(y) # Get variable names for attr
      return(output)
      }

# Summary function for VAR models
"summary.VAR" <- function(object, ...)
{
    labels <- c("Reduced Form VAR", "Constants",
                "Posterior Regression Coefficients",
                "Posterior error covariance")
    rf.list <- list(labels=c("Number of observations","Degrees of freedom per equation"),
                    values=c(nrow(object$Y), nrow(object$Y)-nrow(object$Bhat)))
    constants.list <- list(labels=c("Constants"),
                           values=round(object$intercept,6))
    ar.list <- list(labels=c("Autoregressive Matrices"),
                    values=round(object$ar,6))
    if(is.na(object$exog[1])==FALSE){
        exog.list <- list(labels=c("Exogenous variable posterior coefficients"),
                          values=round(object$exog,6))
    } else { exog.list <- list(labels=NULL,values=NULL) }
    vcv.list <- list(labels=c("Posterior error covariance"),
                     values=object$mean.S)
    values <- list(rf=rf.list, ar=ar.list, constants=constants.list,
                   exog=exog.list, vcv=vcv.list)
    output <- list(labels=labels, values=values)
    list.print(output)
}
