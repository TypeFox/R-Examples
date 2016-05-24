"MleCslogistic"<-function(formula,type=TRUE,intercept=TRUE,method="BFGS",maxiter=1000,data,...)
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
	loglike <- 0
    
	if(type)
	{  
		  dimx <- p+dima
		  coeff <- rep(0,dimx)
        	
        	neg.loglike <- function(coeff)
        	{
        		.Fortran("cloga",as.integer(n),as.integer(nrec),as.integer(p),
        	         as.integer(dima),as.integer(y),as.double(x),
        	         as.double(coeff),as.integer(perm),loglike=as.double(loglike),PACKAGE="cslogistic")$loglike
			}
                
	        a <- optim(par=rep(1.2,dimx),fn=neg.loglike,method=method,hessian = TRUE,control=list(maxit=maxiter))
	        pnames <- c(dimnames(x)[[2]])
	        for(i in 1:(n-1))
		    {
	           for(j in (i+1):n)
			   {
	                bn <- paste("alpha",i,j,sep="")
	           	    pnames <- c(pnames,bn)
	           }
	        }
        }        
                
        else
	    {  
        
			if(intercept)
			{
				dimx <- n+(p-1)+dima
				coeff <- rep(0,dimx)
				neg.loglike <- function(coeff)
				{
					.Fortran("clogc",as.integer(n),as.integer(nrec),as.integer(p),
							 as.integer(dima),as.integer(y),as.double(x),
							 as.double(coeff),as.integer(perm),loglike=as.double(loglike),PACKAGE="cslogistic")$loglike
				}
				a <- optim(rep(1.2,dimx),fn=neg.loglike,method=method,hessian = TRUE,control=list(maxit=maxiter))
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
                
			else
			{
				dimx <- n*p+dima
				coeff <- rep(0,dimx)
				neg.loglike <- function(coeff)
				{
						.Fortran("clogb",as.integer(n),as.integer(nrec),as.integer(p),
								 as.integer(dima),as.integer(y),as.double(x),
								 as.double(coeff),as.integer(perm),loglike=as.double(loglike),PACKAGE="cslogistic")$loglike
				}
                
				a <- optim(rep(1.2,dimx),fn=neg.loglike,method=method,hessian = TRUE,control=list(maxit=maxiter))
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

		coeff <- a$par
		names(coeff) <- pnames
	
        loglike <- -a$value  
        
        conv <- -a$convergence
        covb <- -solve(-a$hessian)
        colnames(covb) <- pnames
        rownames(covb) <- pnames

        se <- sqrt(diag(covb)) 
        names(se) <- pnames
        
        corr <- covb/outer(se, se)
        colnames(corr) <- pnames
        rownames(corr) <- pnames
        
        tvalue <- coeff/se
        
        model.name <- "Conditional specified logistic regression model"
        z <- list(modelname = model.name,call=cl,coefficients = coeff,se=se,tvalue=tvalue,nvar=n,method=method,
        	  loglike=loglike,np=dimx,cov=covb,corr=corr,nrec=nrec)
        class(z) <- c("MleCslogistic")
        
        if(conv!=0)stop("Convergence criteria not achieved, code # ",conv)
        z
}

