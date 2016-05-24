
"BayesCslogistic" <- function(formula, type = TRUE, intercept = TRUE, burnin = 1000, mcmc = 10000,
           thin=1, tune=1.1, beta.start = NA, b0 = 0, B0 = 0, ...)
{
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data","na.action"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	
	y <- model.response(mf,"numeric")
	n <- ncol(y)
	nrec <- nrow(y)
	x <- as.matrix(model.matrix(formula))
	p <- ncol(x)
	dima <- n*(n-1)/2
	perm <- matrix(0,nrow=(2**n),ncol=n)
	arate <- 0
	nstore <- mcmc/thin
	sam <- sample(1:29000,3)
        
 ## seeds
	seed1 <- sam[1]
	seed2 <- sam[2]
    seed3 <- sam[3]
	
 ## equal effect model
    if(type)
	{  
      ## dimension
                
         dimen <- p+dima
                
	  ## maximum likelihood estimates
                
         fit0 <- MleCslogistic(formula)
                
	  ## starting values
                
         if(is.na(beta.start)) beta.start <- fit0$coeff

      ## prior
                
         if(b0==0) b0 <- rep(0,dimen)
		 if(B0==0) B0 <- diag(rep(1000,dimen))
    
      ## covariance matrix for metropolis algorithm
                
         propv <- diag(tune,dimen)%*%(solve(solve(fit0$cov)+solve(B0)))%*%diag(tune,dimen)

	  ## working variables
                
		 betac <- rep(0,dimen)
         workv1 <- rep(0,dimen)
         workv2 <- rep(0,dimen)
         workv3 <- rep(0,dimen)
         workm1 <- rep(0,(dimen*(dimen+1)/2))
         workm2 <- matrix(0,nrow=dimen,ncol=dimen)
         workm3 <- matrix(0,nrow=dimen,ncol=dimen)
         iflag <- rep(0,dimen)
         storemat <- matrix(0,nrow=nstore,ncol=dimen)

	  ## call the fortran function to fit the model

         a <- .Fortran("mha",mcmc=as.integer(mcmc),burnin=as.integer(burnin),thin=as.integer(thin),
                        nstore=as.integer(nstore),as.integer(seed1),as.integer(seed2),as.integer(seed3),
                        propv=as.double(propv),p=as.integer(p),dima=as.integer(dima),dimen=as.integer(dimen),
						n=as.integer(n),nrec=as.integer(nrec),y=as.integer(y),x=as.double(x),as.integer(perm),
                        beta.start=as.double(beta.start),as.double(betac),b0=as.double(b0),
						B0=as.double(B0),mat=as.double(storemat),arate=as.double(arate),
                        as.double(workm1),as.double(workm2),as.double(workm3),
                        as.double(workv1),as.double(workv2),as.double(workv3),
                        as.integer(iflag),PACKAGE="cslogistic")

	  ## save names of original covariates
                
	     pnames <- c(dimnames(x)[[2]])
		 for(i in 1:(n-1))
		 {
			 for(j in (i+1):n)
			 {
			     bn <- paste("alpha",i,j,sep="")
	           	 pnames<-c(pnames,bn)
			 }
		  }
	}        

 ## different effect models               
    else
	{  
        
	  ## different intercept only
                
         if(intercept)
		 {
		   ## dimension  
                     
			  dimen <- n+(p-1)+dima

		   ## maximum likelihood estimates
                
              fit1 <- MleCslogistic(formula,type=type)
     
		   ## starting values
                  
              if(is.na(beta.start)) beta.start <- fit1$coeff

		   ## prior
                
              if(b0==0) b0 <- rep(0,dimen)
              if(B0==0) B0 <- diag(rep(1000,dimen))
       
		   ## covariance matrix for metropolis algorithm
                
              propv <- diag(tune,dimen)%*%(solve(solve(fit1$cov)+solve(B0)))%*%diag(tune,dimen)

		   ## working variables
                
              betac <- rep(0,dimen)
              workv1 <- rep(0,dimen)
              workv2 <- rep(0,dimen)
              workv3 <- rep(0,dimen)
              workm1 <- rep(0,(dimen*(dimen+1)/2))
              workm2 <- matrix(0,nrow=dimen,ncol=dimen)
              workm3 <- matrix(0,nrow=dimen,ncol=dimen)
              iflag <- rep(0,dimen)
              storemat <- matrix(0,nrow=nstore,ncol=dimen)

		   ## call the fortran function to fit the model

              a <- .Fortran("mhc",mcmc=as.integer(mcmc),burnin=as.integer(burnin),thin=as.integer(thin),
                                       nstore=as.integer(nstore),as.integer(seed1),as.integer(seed2),as.integer(seed3),
                                       propv=as.double(propv),p=as.integer(p),dima=as.integer(dima),dimen=as.integer(dimen),
                                       n=as.integer(n),nrec=as.integer(nrec),y=as.integer(y),x=as.double(x),as.integer(perm),
                                       beta.start=as.double(beta.start),as.double(betac),b0=as.double(b0),
                                       B0=as.double(B0),mat=as.double(storemat),arate=as.double(arate),
                                       as.double(workm1),as.double(workm2),as.double(workm3),
                                       as.double(workv1),as.double(workv2),as.double(workv3),
                                       as.integer(iflag),PACKAGE="cslogistic")
  
		   ## save names of original covariates       
     
              pnames <- c(dimnames(x)[[2]])                     
              pnames0 <- paste(pnames[1],seq(1:n),sep=":")
              pnames1 <- pnames[-1]
              pnames <- c(pnames0,pnames1)
                     
			  for(i in 1:(n-1))
			  {
		          for(j in (i+1):n)
				  {
					  bn <- paste("alpha",i,j,sep="")
				      pnames <- c(pnames,bn)
				  }
			   }
		}

	 ## different effects in all the covariates
                
        else
		{
		   ## dimension  
                     
              dimen <- n*p+dima

		   ## maximum likelihood estimates
                
              fit2 <- MleCslogistic(formula,type=type,intercept=intercept)
                     
		   ## starting values
                  
              if(is.na(beta.start)) beta.start <- fit2$coeff

		   ## prior
                
              if(b0==0) b0 <- rep(0,dimen)
              if(B0==0) B0 <- diag(rep(1000,dimen))
       
		   ## covariance matrix for metropolis algorithm
                
              propv <- diag(tune,dimen)%*%(solve(solve(fit2$cov)+solve(B0)))%*%diag(tune,dimen)

		   ## working variables
                
              betac <- rep(0,dimen)
              workv1 <- rep(0,dimen)
              workv2 <- rep(0,dimen)
              workv3 <- rep(0,dimen)
              workm1 <- rep(0,(dimen*(dimen+1)/2))
              workm2 <- matrix(0,nrow=dimen,ncol=dimen)
              workm3 <- matrix(0,nrow=dimen,ncol=dimen)
              iflag <- rep(0,dimen)
              storemat <- matrix(0,nrow=nstore,ncol=dimen)

		   ## call the fortran function to fit the model

			  a <- .Fortran("mhb",mcmc=as.integer(mcmc),burnin=as.integer(burnin),thin=as.integer(thin),
                                       nstore=as.integer(nstore),as.integer(seed1),as.integer(seed2),as.integer(seed3),
                                       propv=as.double(propv),p=as.integer(p),dima=as.integer(dima),dimen=as.integer(dimen),
                                       n=as.integer(n),nrec=as.integer(nrec),y=as.integer(y),x=as.double(x),as.integer(perm),
                                       beta.start=as.double(beta.start),as.double(betac),b0=as.double(b0),
                                       B0=as.double(B0),mat=as.double(storemat),arate=as.double(arate),
                                       as.double(workm1),as.double(workm2),as.double(workm3),
                                       as.double(workv1),as.double(workv2),as.double(workv3),
                                       as.integer(iflag),PACKAGE="cslogistic")
  
		   ## save names of original covariates       
                     
              pnames <- c(dimnames(x)[[2]])
              pnames <- paste(pnames,1,sep=":")
                     
			  for(i in 2:n)
			  {
				  nbv <- paste(dimnames(x)[[2]],":",i,sep="")
				  pnames <- c(pnames,nbv)
			  }
              for(i in 1:(n-1))
			  {
				  for(j in (i+1):n)
				  {
					  bn <- paste("alpha",i,j,sep="")
					  pnames <- c(pnames,bn)
				  }
			   }
			}      
        }        

        mat <- matrix(a$mat,nstore,dimen,byrow=FALSE)
        colnames(mat) <- pnames
        
        coeff <- rep(0,dimen)
        for(i in 1:dimen)
	    {
            coeff[i] <- mean(mat[,i])
        }
        names(coeff) <- pnames

        arate <- a$arate
 
        
        model.name <- "Bayesian conditionally specified logistic regression model"
        z <- list(modelname = model.name,call=cl,coefficients = coeff,mcmc=mcmc,burnin=burnin,
                  thin=thin,mat=mat,arate=arate,dimen=dimen,pnames=pnames)
        class(z) <- c("BayesCslogistic")
        z
}



