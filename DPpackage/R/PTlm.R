### PTlm.R                   
### Fit a semiparametric regression model.
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 13-06-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author's contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###

"PTlm"<-
function(formula,ngrid=200,grid=NULL,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PTlm")

"PTlm.default"<-
function(formula,
         ngrid=200,
         grid=NULL,
         prior,
         mcmc,
         state,
         status, 
         data=sys.frame(sys.parent()),
         na.action=na.fail)
{
         #########################################################################################
         # call parameters
         #########################################################################################

	       cl <- match.call()
	       mf <- match.call(expand.dots = FALSE)
	       m <- match(c("formula", "data","na.action"), names(mf), 0)
	       mf <- mf[c(1, m)]
	       mf$drop.unused.levels <- TRUE
	       mf[[1]] <- as.name("model.frame")
	       mf <- eval(mf, parent.frame())

         #########################################################################################
         # data and model structure
         #########################################################################################
           if(is.null(prior$frstlprob))
           {
              mdzero <- 1
		   }
           else
           {
              mdzero <- 0
              if(prior$frstlprob)mdzero <- 1
		   }

     	   y <- model.response(mf,"numeric")
  	       nrec <- length(y)
		   x <-  model.matrix(formula)
  	       p <- ncol(x)

		   if(mdzero==0)
		   {
              p <- p-1
			  if(p>1)
			  {
                 x <- x[,-1]
              }   
           }   
  	       nfixed<-p
           if(p==0)
           {
              nfixed <- 0
              p <- 1
              x <- matrix(0,nrow=nrec,ncol=1)
           }

         #########################################################################################
         # Elements for Pseudo Countour Probabilities' computation
         #########################################################################################
           Terms <- if (missing(data)) 
                terms(formula)
           else terms(formula, data = data)

           ntanova<-0   
           if((mdzero==1) & (nfixed==1)) ntanova<-1

           possiP <- NULL
           if((nfixed>0) & (ntanova==0))
           {
              mat <- attr(Terms,"factors")
              namfact <- colnames(mat)
              nvar <- dim(mat)[1]
              nfact <- dim(mat)[2]
              possiP <- matrix(0,ncol=2,nrow=nfact)
              if (missing(data)) dataF <- model.frame(formula=formula,xlev=NULL)
               dataF <- model.frame(formula=formula,data,xlev=NULL)
              namD <- names(dataF)
              isF <- sapply(dataF, function(x) is.factor(x) || is.logical(x))
              nlevel <- rep(0,nvar)
              for(i in 1:nvar)
              {
                if(isF[i])
                {
                   nlevel[i]<-length(table(dataF[[i]]))
                }
                else
                {
                   nlevel[i]<-1
                }
              }

              startp<-1+1
              if(mdzero==1)startp<-1+2

              for(i in 1:nfact)
              {
                tmp1<-1
                for(j in 1:nvar)
                {
                    if(mat[j,i]==1 && isF[j])
                    {
                       tmp1<-tmp1*(nlevel[j]-1)
                    }
                }
                endp<-startp+tmp1-1
                possiP[i,1]<-startp    
                possiP[i,2]<-endp
                startp<-endp+1
              }
              dimnames(possiP)<-list(namfact,c("Start","End"))
           }   

         #########################################################################################
         # prior information
         #########################################################################################

           if(nfixed==0)
           {
              betapv <- matrix(0,nrow=1,ncol=1)
              betapm <- rep(0,1)
              propv <- matrix(0,nrow=1,ncol=1)
           }
           else
           {
              betapm <- prior$beta0
              betapv <- prior$Sbeta0
              propv <- solve(t(x)%*%x+solve(prior$Sbeta0))

              if(length(betapm)!=p)
              { 
                     stop("Error in the dimension of the mean of the normal prior for the fixed effects.\n")     
              }

              if(dim(betapv)[1]!=p || dim(betapv)[2]!=p)
              { 
                     stop("Error in the dimension of the covariance of the normal prior for the fixed effects.\n")     
              }

           }

  	       if(is.null(prior$a0))
		   {
			  alpharand<-0 
			  a0b0<-c(-1,-1)
			  alpha<-prior$alpha
		   }
		   else
		   {
			   alpharand<-1 	 
			   a0b0<-c(prior$a0,prior$b0)
			   alpha<-rgamma(1,prior$a0,prior$b0)
		   }
	
		   if(mdzero==1)
           {
			  murand<-0
              mu<-0
              m0<-0
              s0<--1
           }
           else
           {
            if(is.null(prior$mub))
            {
               murand<-0 
               if(is.null(prior$mu))
               { 
                  stop("*mu* must be specified in the prior object when it is not considered as random.\n")     
               }
               if(length(prior$mu) != 1)
               { 
                  stop("Error in the dimension of the mean the centering distribution.\n")     
               }
               m0<-1
               s0<--1
            }
            else
            {
               murand<-1
               s0<-prior$Sb
               m0<-prior$mub
               if(length(prior$mub) != 1)
               { 
                  stop("Error in the dimension of the mean of the normal prior for the mean of the centering distribution.\n")     
               }
               if(!is.null(dim(s0)) && ( dim(s0)[1]!=1 || dim(s0)[2]!=1 ))
               { 
                  stop("Error in the dimension of the covariance of the normal prior for the mean of the centering distribution.\n")     
               }
             }
           }


      	   if(is.null(prior$tau1))
  	       {
  	          sigmarand <- 0
              if(is.null(prior$sigma2))
              { 
                 stop("The variance *sigma2* must be specified in the prior object when it is not considered as random.\n")     
              }
              sigma2 <- prior$sigma2
              tau1<- -1
              tau2<- -1
		   }
           else
           {
              sigmarand <- 1
              tau1 <- prior$tau1
              tau2 <- prior$tau2

              if(tau1<=0 || tau2<0)
              { 
                 stop("The hyperparameter of the gamma prior for the centering variance must be positive")     
              }
              sigma2 <- tau1/tau2
           }

 	       tau <- c(tau1,tau2)

           maxm <- prior$M
	
         #########################################################################################
         # mcmc specification
         #########################################################################################

           if(missing(mcmc))
           {
              nburn <- 1000
              nsave <- 1000
              nskip <- 0
              ndisplay <- 100
              mcmcvec <- c(nburn,nskip,ndisplay)
           }
           else
           {
              mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
              nsave <- mcmc$nsave
           }

         #########################################################################################
         # output
         #########################################################################################
         
           acrate<-rep(0,4)
           if(mdzero==1)
           {
  	          fit0  <- glm.fit(x, y, family= gaussian(link = "identity"))   
  	          beta <- coefficients(fit0)
			  left <- min(resid(fit0))-0.5*sqrt(var(resid(fit0)))
              right <- max(resid(fit0))+0.5*sqrt(var(resid(fit0)))
		   }
           else
           {
  	          fit0 <- glm.fit(cbind(x,rep(1,nrec)), y, family= gaussian(link = "identity"))   
  	          beta <- coefficients(fit0)[1:p]
              left <- min(coefficients(fit0)[p+1]+resid(fit0))-0.5*sqrt(var(resid(fit0)))
              right <- max(coefficients(fit0)[p+1]+resid(fit0))+0.5*sqrt(var(resid(fit0)))
           }
 
           cpo <- rep(0,nrec)
         
           if(is.null(grid))
	       {
               grid <- seq(left,right,length=ngrid)
           }
           else
           {
               ngrid <- length(grid)  
		   }

   	       f <- rep(0,ngrid)
		   thetasave <- matrix(0, nrow=nsave, ncol=nfixed+3)
		   randsave  <- matrix(0, nrow=nsave, ncol=nrec+1)

         #########################################################################################
         # parameters depending on status
         #########################################################################################

   	       if(status==TRUE)
	       {
	          if(alpharand==1)
	          {
	             alpha <- 1
	          }
	          else
	          {
	             alpha <- prior$alpha  
	          }

  	          if(murand==1)
	          {
				 mu <- coefficients(fit0)[p+1]
  	          }
  	          else
  	          {
				if(mdzero==0)mu <- prior$mu
  	          }
  	        
  	          if(sigmarand==1)
	          {
                   sigma2 <- sum((resid(fit0))**2)/(nrec-p)
			  }
			  else
			  {
                   sigma2 <- prior$sigma2
			  }

			  if(mdzero==1)
			  {
				 v <- resid(fit0)
       		  }
  		      else
  		      {
                   v <- mu+resid(fit0)
			  }   
	        
			  mcmcad <- c(0.1,0.0,
                          -6.0,1.0,0.0,0.0,0.0,
                          -6.0,1.0,0.0,0.0,0.0,
                          -3.0,1.0,0.0,0.0,0.0,
                          0.0)
  		   }	

           if(status==FALSE)
	       {
			   beta <- state$beta
	        
	           if(alpharand==1)
	           {
	              alpha <- state$alpha 
	           }
	           if(murand==1)
	           {
	              mu <- state$mu
	           }
	           if(sigmarand==1)
	           {
	              sigma2 <- state$sigma2
	           }

	           v <- state$v
	        
	           if(is.null(state$mcmcad))
	           {
                   mcmcad <- c(0.1,0.0,
                            -6.0,1.0,0.0,0.0,0.0,
                            -6.0,1.0,0.0,0.0,0.0,
                            -3.0,1.0,0.0,0.0,0.0,                            
                             0.0)
	           }
	           else
	           {
	              mcmcad <- state$mcmcad
	           }
	       }

         #########################################################################################
         # working space
         #########################################################################################
           seed <- c(sample(1:29000,1),sample(1:29000,1))
           mdzero <- 0
 	       iflag <- rep(0,p)
 	       whicho <- rep(0,nrec)
	       whichn <- rep(0,nrec)
	       betac <- rep(0,p)
		   workm1 <- matrix(0,nrow=p,ncol=p)
	       workm2 <- matrix(0,nrow=p,ncol=p)
		   workmh1 <- rep(0,(p*(p+1)/2))
	       workv1 <- rep(0,p)
	       workv2 <- rep(0,p)
 	       vc <- rep(0,nrec)
	       xtx <- matrix(0,nrow=p,ncol=p)
	
         #########################################################################################
         # calling the fortran code
         #########################################################################################

           if(is.null(prior$M))
           {
               foo <- .Fortran("ptlm",
					mdzero    =as.integer(mdzero),
					ngrid     =as.integer(ngrid),
					nrec      =as.integer(nrec),
					p         =as.integer(p),
					x         =as.double(x),	 	
					y         =as.double(y),
					a0b0      =as.double(a0b0),
					betapm    =as.double(betapm),		
					betapv    =as.double(betapv),		
					tau       =as.double(tau),
					m0        =as.double(m0),
					s0        =as.double(s0), 
					mcmc      =as.integer(mcmcvec),
					nsave     =as.integer(nsave),
					propv     =as.double(propv),
					mcmcad    =as.double(mcmcad),
					seed      =as.integer(seed),
					acrate    =as.double(acrate),
					randsave  =as.double(randsave),
					thetasave =as.double(thetasave),
					cpo       =as.double(cpo),
					f         =as.double(f),
					alpha     =as.double(alpha),		
					beta      =as.double(beta),
					mu        =as.double(mu),
					sigma2    =as.double(sigma2),
					v         =as.double(v),
					betac     =as.double(betac),
					iflag     =as.integer(iflag),
					vc        =as.double(vc),
					workm1    =as.double(workm1),
					workm2    =as.double(workm2),
					workmh1   =as.double(workmh1),
					workv1    =as.double(workv1),
					workv2    =as.double(workv2),
					grid      =as.double(grid),
					whicho    =as.integer(whicho),
					whichn    =as.integer(whichn),
					xtx       =as.double(xtx),
					PACKAGE="DPpackage")	                  
         }
         else
         {
            foo <- .Fortran("ptlmp",
					maxm      =as.integer(maxm),
					mdzero    =as.integer(mdzero),
					ngrid     =as.integer(ngrid),
					nrec      =as.integer(nrec),
					p         =as.integer(p),
					x         =as.double(x),	 	
					y         =as.double(y),
					a0b0      =as.double(a0b0),
					betapm    =as.double(betapm),		
					betapv    =as.double(betapv),		
					tau       =as.double(tau),
					m0        =as.double(m0),
					s0        =as.double(s0), 
					mcmc      =as.integer(mcmcvec),
					nsave     =as.integer(nsave),
					propv     =as.double(propv),
					mcmcad    =as.double(mcmcad),
					seed      =as.integer(seed),
					acrate    =as.double(acrate),
					randsave  =as.double(randsave),
					thetasave =as.double(thetasave),
					cpo       =as.double(cpo),
					f         =as.double(f),
					alpha     =as.double(alpha),		
					beta      =as.double(beta),
					mu        =as.double(mu),
					sigma2    =as.double(sigma2),
					v         =as.double(v),
					betac     =as.double(betac),
					iflag     =as.integer(iflag),
					vc        =as.double(vc),
					workm1    =as.double(workm1),
					workm2    =as.double(workm2),
					workmh1   =as.double(workmh1),
					workv1    =as.double(workv1),
					workv2    =as.double(workv2),
					grid      =as.double(grid),
					whicho    =as.integer(whicho),
					whichn    =as.integer(whichn),
					xtx       =as.double(xtx),
					PACKAGE="DPpackage")	         
           }
         
         #########################################################################################
         # save state
         #########################################################################################
	
 	       thetasave <- matrix(foo$thetasave,nrow=nsave, ncol=(nfixed+3))
 	       randsave <- matrix(foo$randsave,nrow=nsave, ncol=(nrec+1))
 	 
  	       colnames(thetasave) <- c(dimnames(x)[[2]],"mu","sigma2","alpha")

           qnames <- NULL
           for(i in 1:nrec)
		   {
             idname <- paste("(Subject",i,sep="=")
             idname <- paste(idname,")",sep="")
             qnames <- c(qnames,idname)
           }
           qnames <- c(qnames,"Prediction")
         
           colnames(randsave) <- qnames

  	       model.name <- "Bayesian Semiparametric Regression Model"
           coeff <- apply(thetasave,2,mean)		
	
	       state <- list(alpha=foo$alpha,
                         beta=foo$beta,
                         v=foo$v,
                         mu=foo$mu,
                         sigma2=foo$sigma2,
     	                 mcmcad=foo$mcmcad)	
		      
	       save.state <- list(thetasave=thetasave,randsave=randsave)

           if(is.null(prior$a0))
	       {
              acrate <- foo$acrate[1:3]
		   }
           else
           {
              acrate <- foo$acrate
           }   

	       z <- list(modelname=model.name,
					 coefficients=coeff,
					 acrate=acrate,
					 call=cl,
				     prior=prior,
					 mcmc=mcmc,
					 state=state,
					 save.state=save.state,
					 cpo=foo$cpo,
				     nrec=nrec,
					 p=p,
					 dens=foo$f,
					 grid=grid,
					 possiP=possiP)
	
		  cat("\n\n")
	      class(z) <- c("PTlm")
	      z 
}


###
### Tools for PTlm: anova, print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2007
### Last modification: 14-08-2007.


"anova.PTlm"<-function(object, ...)
{

######################################################################################
cregion<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a simultaneous credible region for a vector 
#  parameter from the MCMC sample
# 
#  Reference: Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995)
#             Bayesian computation and stochastic systems (with Discussion)
#             Statistical Science, vol. 10, 3 - 66, page 30
#  and        Held, L. (2004) Simultaneous inference in risk assessment; a Bayesian 
#             perspective In: COMPSTAT 2004, Proceedings in Computational 
#             Statistics (J. Antoch, Ed.) 213 - 222, page 214
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are computed.
######################################################################################
{
    #Basic information
     nmonte<-dim(x)[1]
     p<-dim(x)[2]
     
    #Ranks for each component
     ranks <- apply(x, 2, rank, ties.method="first")
     
    #Compute the set S={max(nmonte+1-min r_i(t) , max r_i(t)): t=1,..,nmonte}
     left <- nmonte + 1 - apply(ranks, 1, min)
     right <- apply(ranks, 1, max)
     S <- apply(cbind(left, right), 1, max)
     S <- S[order(S)]
    
    #Compute the credible region
     k <- floor(nmonte*probs)     
     tstar <- S[k]
     out<-list()
     for(i in 1:length(tstar))
     {
        upelim <- x[ranks == tstar[i]]
        lowlim <- x[ranks == nmonte + 1 - tstar[i]]    
        out[[i]] <- rbind(lowlim, upelim)
        rownames(out[[i]]) <- c("Lower", "Upper")
        colnames(out[[i]]) <- colnames(x)
     }
     names(out) <- paste(probs)
     return(out)
}

######################################################################################
cint<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a credible interval from the MCMC sample
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are to be computed.
######################################################################################
{
    #Compute the credible interval
     delta<-(1-probs)/2
     lprobs<-cbind(delta,probs+delta) 
     out<-matrix(quantile(x,probs=lprobs),ncol=2)
     colnames(out) <- c("Lower","Upper")
     rownames(out) <- paste(probs)
     return(out)
}

######################################################################################
hnulleval<-function(mat,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     npar<-dim(mat)[2]   
     lower<-rep(0,npar)
     upper<-rep(0,npar)
     for(i in 1:npar)
     {
        lower[i]<-mat[1,i]< hnull[i]
        upper[i]<-mat[2,i]> hnull[i]
     }
     total<-lower+upper
     out<-(sum(total==2) == npar)
     return(out)
}

######################################################################################
hnulleval2<-function(vec,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     lower<-vec[1]< hnull
     upper<-vec[2]> hnull

     total<-lower+upper
     out<-(total==2)
     return(out)
}


######################################################################################
pcp<-function(x,hnull=NULL,precision=0.001,prob=0.95,digits=digits)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Region)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-rep(0,dim(x)[2])
    if (dim(x)[2]!=length(hnull)) stop("Dimension of x and hnull must be equal!!")

    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cregion(x,probs=probsf)

    is.hnull <- hnulleval(cr[[2]],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval(cr[[length(cr)]],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval(cr[[i]],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1 - probs[first]
       }
    }
    output <- list(cr=cr[[1]], prob=prob, pval=pval,hnull=hnull)
    return(output)
}


######################################################################################
pcp2<-function(x,hnull=NULL,precision=0.001,prob=0.95)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Interval)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-0
    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cint(x,probs=probsf)

    is.hnull <- hnulleval2(cr[2,],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval2(cr[(neval+1),],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval2(cr[i,],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1-probs[first]
       }
    }
    output <- list(cr=cr[1,], prob=prob, pval=pval,hnull=hnull)
    return(output)
}

######################################################################################
######################################################################################
######################################################################################

    possiP<-object$possiP
    nfact<-dim(possiP)[1]
    P<-rep(0,nfact)
    df<-rep(0,nfact)
    
    for(i in 1:nfact)
    {
        df[i]<-1
        if((possiP[i,2]-possiP[i,1])>0)
        { 
           x<-matrix(object$save.state$thetasave[,possiP[i,1]:possiP[i,2]])
           foo<-pcp(x=x) 
           P[i]<-foo$pval
           df[i]<-(possiP[i,2]-possiP[i,1])+1
        }
        else
        {
           x<-object$save.state$thetasave[,possiP[i,1]:possiP[i,2]]
           foo<-pcp2(x=x) 
           P[i]<-foo$pval
        }
    }

    table <- data.frame(df,P) 
    dimnames(table) <- list(rownames(possiP), c("Df","PsCP"))
    structure(table, heading = c("Table of Pseudo Contour Probabilities\n", 
        paste("Response:", deparse(formula(object)[[2]]))), class = c("anovaPsCP",
        "data.frame"))
}




"print.PTlm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    cat("\nPosterior Inference of Parameters:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"plot.PTlm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
{

fancydensplot<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
# Author: AJV, 2006
#
{
	dens <- density(x,n=npts)
	densx <- dens$x
	densy <- dens$y

	meanvar <- mean(x)
	densx1 <- max(densx[densx<=meanvar])
	densx2 <- min(densx[densx>=meanvar])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ymean <- densy1 + ((densy2-densy1)/(densx2-densx1))*(meanvar-densx1)
        

        if(hpd==TRUE)
	{
		alpha<-0.05
		alow<-rep(0,2)
        	aupp<-rep(0,2)
        	n<-length(x)
		a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(x),
		                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
		xlinf<-a$alow[1]            
		xlsup<-a$aupp[1]            
	}
	else
	{
		xlinf <- quantile(x,0.025)
		xlsup <- quantile(x,0.975)
	}

	densx1 <- max(densx[densx<=xlinf])
	densx2 <- min(densx[densx>=xlinf])
	densy1 <- densy[densx==densx1]
	densy2 <- densy[densx==densx2]
	ylinf <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlinf-densx1)

	densx1 <- max(densx[densx<=xlsup])
	densx2 <- min(densx[densx>=xlsup])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ylsup <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlsup-densx1)

        plot(0.,0.,xlim = c(min(densx), max(densx)), ylim = c(min(densy), max(densy)),
             axes = F,type = "n" , xlab=xlab, ylab=ylab, main=main, cex=1.2)

        
        xpol<-c(xlinf,xlinf,densx[densx>=xlinf & densx <=xlsup],xlsup,xlsup)
        ypol<-c(0,ylinf,densy[densx>=xlinf & densx <=xlsup] ,ylsup,0)
             
        polygon(xpol, ypol, border = FALSE,col=col)
        
        lines(c(min(densx), max(densx)),c(0,0),lwd=1.2)
        
        segments(min(densx),0, min(densx),max(densy),lwd=1.2)
        
        lines(densx,densy,lwd=1.2)
             
        segments(meanvar, 0, meanvar, ymean,lwd=1.2)
        segments(xlinf, 0, xlinf, ylinf,lwd=1.2)
        segments(xlsup, 0, xlsup, ylsup,lwd=1.2)

	axis(1., at = round(c(xlinf, meanvar,xlsup), 2.), labels = T,pos = 0.)
        axis(1., at = round(seq(min(densx),max(densx),length=15), 2.), labels = F,pos = 0.)
        axis(2., at = round(seq(0,max(densy),length=5), 2.), labels = T,pos =min(densx))
}


    if(is(x, "PTlm")){
        if(is.null(param))
	{
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:(n-1)){
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
	       {
	          hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
	       }
	       else
	       {
                  fancydensplot(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }   
           }
           
           if(is.null(x$prior$a0))
           {
              cat("")
           }
           else
           {
               title1<-paste("Trace of",pnames[n],sep=" ")
               title2<-paste("Density of",pnames[n],sep=" ")       
               plot(x$save.state$thetasave[,n],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-c("Predictive Error Density")
           plot(x$grid,x$dens,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="values")
           
        }
        else
        {
            coef.p<-x$coefficients
	    n<-length(coef.p)
	    pnames<-names(coef.p)
	    poss<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }

            if (poss==0 && param !="predictive") 
	    {
	      stop("This parameter is not present in the original model.\n")
	    }

            if (param !="predictive") 
            {
	       par(ask = ask)
	       layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
               title1<-paste("Trace of",pnames[poss],sep=" ")
               title2<-paste("Density of",pnames[poss],sep=" ")       
               plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            else
            {
               title1<-c("Predictive Error Density")
               plot(x$grid,x$dens,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="values")
            }
        }
   }
}


"summary.PTlm"<-function(object, hpd=TRUE, ...) 
{
    dimen<-object$p+2
    coef.p<-object$coefficients[1:dimen]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])
    
    alpha<-0.05
    
    for(i in 1:dimen){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,i]))
        coef.m[i]<-median(object$save.state$thetasave[,i])
        vec<-object$save.state$thetasave[,i]
        n<-length(vec)
        
        if(hpd==TRUE)
        {
        
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                coef.l[i]<-a$alow[1]            
                coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    
    if(hpd==TRUE)
    {
       dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
    }            
    else
    {
       dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
    }
    
    
    ans <- c(object[c("call", "modelname")])
    
    ans$coefficients<-coef.table

    ans$cpo<-object$cpo
    

    if(is.null(object$prior$a0))
    {
	
    }
    else
    {
         dimen<-1
         coef.p<-object$coefficients[(object$p+2)]
         coef.sd<-rep(0,dimen)
         coef.se<-rep(0,dimen)
         coef.l<-rep(0,dimen)
         coef.u<-rep(0,dimen)
         coef.m<-rep(0,dimen)
         names(coef.sd)<-names(object$coefficients[(object$p+2)])
         names(coef.l)<-names(object$coefficients[(object$p+2)])
         names(coef.u)<-names(object$coefficients[(object$p+2)])
         for(i in 1:dimen){
             alow<-rep(0,2)
             aupp<-rep(0,2)
             coef.sd[i]<-sqrt(var(object$save.state$thetasave[,(object$p+1+i)]))
             coef.m[i]<-median(object$save.state$thetasave[,(object$p+1+i)])
             vec<-object$save.state$thetasave[,(object$p+i)]
             n<-length(vec)
        
             if(hpd==TRUE)
             {
                 a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                              alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                 coef.l[i]<-a$alow[1]            
                 coef.u[i]<-a$aupp[1]            
             }
             else
             {
                 coef.l[i]<-quantile(vec,0.025) 
                 coef.u[i]<-quantile(vec,0.975) 
             }
         }
      
         coef.se<-coef.sd/sqrt(n)
         coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

         if(hpd==TRUE)
         {
              dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
         }            
         else
         {
              dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
         }

         ans$prec<-coef.table
    }


    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryPTlm"
    return(ans)
}



"print.summaryPTlm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) 
    {
        cat("\nRegression coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$prec)) 
    {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat(" \n")

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


