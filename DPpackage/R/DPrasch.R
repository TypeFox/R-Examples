### DPrasch.R                   
### Fit a Rasch model with a Dirichlet Process prior
### for the random effect distribution
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 04-09-2009.
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


DPrasch <- function(y,prior,mcmc,offset,state,status,grid=seq(-10,10,length=1000),
                         data=sys.frame(sys.parent()),compute.band=FALSE)
UseMethod("DPrasch")

DPrasch.default<-
function(y,
         prior,
         mcmc,
         offset=NULL,
         state,
         status,
         grid=seq(-10,10,length=1000),
         data=sys.frame(sys.parent()),
         compute.band=FALSE)
{
         #########################################################################################
         # call parameters
         #########################################################################################
           cl <- match.call()
	       y <- as.matrix(y)
	  
         #########################################################################################
         # data structure
         #########################################################################################
     	   nsubject <- nrow(y)
	       p <- ncol(y)

           ywork <- y
          
           datastrm <- NULL
           nmissing <- 0
           total <- 0
          
           for(i in 1:nsubject)
           {
              for(j in 1:p)
              {
                  if(is.na(y[i,j]))
                  {
                     nmissing <- nmissing+1
                     datastrm <- rbind(datastrm,c(i,j))   
                  }
                  else
                  {
                     total <- total+y[i,j]            
                  }
                  
              }
           }
          
           nrec <- nsubject*p-nmissing
          
           if(nmissing>0)
           {
              imiss <- 1 
              for(i in 1:nmissing)
              {
                   ywork[datastrm[i,1],datastrm[i,2]] <- rpois(1,total/nrec)               
              }
           }
           else
           {
              imiss <- 0
              nmissing <- 1
              datastrm <- matrix(0,nrow=1,ncol=2)
           }


         #########################################################################################
         # prior information
         #########################################################################################

       	   if(is.null(prior$a0))
  	       {
  	          a0 <- -1
  	          b0 <- -1 
  	          alpha <- prior$alpha
		  	  alpharand <- 0
  	       }
           else
           {
              a0 <- prior$a0
  	          b0 <- prior$b0
  	          alpha <- 1
			  alpharand <- 1
		   }
  	       a0b0 <- c(a0,b0)
  	 
     	   if(is.null(prior$tau1))
  	       {
              tau1 <- -1
              tau2 <- -1
              sigma2 <- prior$sigma2
              sigmainv <- 1/sigma2
              sigmarand <- 0
  	       }
  	       else
  	       {
              tau1 <- prior$tau1
              tau2 <- prior$tau2
              sigma2 <-1
              sigmainv <- 1/sigma2
  	          sigmarand <- 1
  	       }
  	 
     	   if(is.null(prior$mub))
  	       {
			  psiinv <- -1
  	          smu <- 0
  	          mu <- prior$mu
  	          murand <- 0
  	       }
  	       else
     	   {
			  psiinv <- (1/prior$Sb)
  	          smu <- psiinv*prior$mub
  	          mu <- rnorm(1,prior$mub,sqrt(prior$Sb))
			  murand <- 1
           }     

           b0 <- prior$beta0
           prec <- solve(prior$Sbeta0)
           sb <- prec%*%prior$beta0

           if(dim(prec)[1] != (p-1)) stop("the dimension of beta0 and Sbeta0 must be p-1")

         #########################################################################################
         # mcmc specification
         #########################################################################################
           mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
           nsave <- mcmc$nsave

         #########################################################################################
         # output
         #########################################################################################
           acrate <- rep(0,2)
		   cpo <- matrix(0,nrow=nsubject,ncol=p)
		   cpov <- rep(0,nsubject)
           ngrid <- length(grid)
           cdfsave <- matrix(0,nrow=nsave,ncol=ngrid)
           thetasave <- matrix(0,nrow=nsave,ncol=p+5)
           randsave <- matrix(0,nrow=nsave,ncol=nsubject+1)
         
         #########################################################################################
         # MLE estimation
         #########################################################################################
         
           RaschMLE <- function(y,nitem,nsubject,offset)
           {
             ywork2 <- NULL
             roffset <- NULL
             id <- NULL
             x <- NULL
             count <- 0
             for(i in 1:nsubject)
             {
				 ywork2 <- c(ywork2,y[i,])
				 roffset <- c(roffset,offset[i,])
				 id <- c(id,rep(i,nitem))

                 aa <- diag(-1,nitem)
				 aa[,1] <- 1
                 x <- rbind(x,aa)
             }
             out <- NULL
             library(nlme)
			 library(MASS)
             fit0 <- glmmPQL(ywork2~x-1+offset(roffset),random = ~ 1 | id,family=binomial(logit), verbose = FALSE) 

             beta <- as.vector(fit0$coeff$fixed[2:nitem])
             b <- as.vector(fit0$coeff$fixed[1]+fit0$coeff$random$id)
			 out$beta <- beta
			 out$b <- b
             out$mu <- fit0$coeff$fixed[1]
             out$sigma2 <- getVarCov(fit0)[1]
             return(out)
		   }
	 
		   if(is.null(offset))
	       {
			  roffset <- matrix(0,nrow=nsubject,ncol=p)
		   }
  	       else
		   {
			  roffset <- offset
		   }
	 
		   fit0 <- RaschMLE(ywork,p,nsubject,roffset)


         #########################################################################################
         # parameters depending on status
         #########################################################################################
       
      	   if(status==TRUE)
	       {
                beta <- fit0$beta
                b <- fit0$b
                bclus <- c(b,rep(0,100))
                ncluster <- nsubject
                ss <- seq(1,nsubject)
                if(murand==1) mu <- fit0$mu
				if(sigmarand==1) sigma2 <- fit0$sigma2
		   }
	 
      	   if(status==FALSE)
	       {
	            alpha <- state$alpha
	            b <- state$b
                bclus <- state$bclus 
                beta <- state$beta
                mu <- state$mu
	            ncluster <- state$ncluster
	            sigma2 <- state$sigma2
	            sigmainv <- 1/sigma2
	            ss <- state$ss
	       }    


         #########################################################################################
         # working space
         #########################################################################################
           betac <- rep(0,(p-1))
           ccluster <- rep(0,nsubject) 
           fsavet <- rep(0,ngrid)         
           iflagp <- rep(0,(p-1))
           prob <- rep(0,nsubject+100)
           seed1 <- sample(1:29000,1)
           seed2 <- sample(1:29000,1)
           seed <- c(seed1,seed2)
           workmhp <- rep(0,(p-1)*p/2)
           workvp <- rep(0,p-1)
           xtx <- matrix(0,nrow=p-1,ncol=p-1)
           xty <- rep(0,p-1)
           cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)
           workcpo <- matrix(0,nrow=nsubject,ncol=p)
      
           if(is.null(offset))
		   {
	          roffset <- matrix(0,nrow=nsubject,ncol=p)
  	       }
  	       else
  	       {
  	          roffset <- offset
  	       }

         #########################################################################################
         # calling the fortran code
         #########################################################################################

		   foo <- .Fortran("sprasch",
				datastrm   =as.integer(datastrm),
				imiss      =as.integer(imiss),
				ngrid      =as.integer(ngrid),  	 	
				nmissing   =as.integer(nmissing),
				nsubject   =as.integer(nsubject),
				p          =as.integer(p),
				y          =as.integer(ywork),
				roffset    =as.double(roffset),
				a0b0       =as.double(a0b0),
				b0         =as.double(b0),
				prec       =as.double(prec),
				psiinv     =as.double(psiinv),
				sb         =as.double(sb),
				smu        =as.double(smu),
				tau1       =as.double(tau1),
				tau2       =as.double(tau2),
				mcmc       =as.integer(mcmcvec),
				nsave      =as.integer(nsave),
                acrate     =as.double(acrate),
				cpo        =as.double(cpo),
				cpov       =as.double(cpov),
				cdfsave    =as.double(cdfsave), 		
				randsave   =as.double(randsave),
				thetasave  =as.double(thetasave),
				alpha      =as.double(alpha),		
				b          =as.double(b),		
				bclus      =as.double(bclus),		
				beta       =as.double(beta),		
                mu         =as.double(mu),		 		
				ncluster   =as.integer(ncluster),
				sigma2     =as.double(sigma2),
				sigmainv   =as.double(sigmainv),
				ss         =as.integer(ss),
				betac      =as.double(betac),		
				ccluster   =as.integer(ccluster),
                cstrt	   =as.integer(cstrt),
                fsavet     =as.double(fsavet), 		
				iflagp     =as.integer(iflagp),
				prob       =as.double(prob),
				seed       =as.integer(seed),
				workmhp    =as.double(workmhp),
				workvp     =as.double(workvp),
				xtx        =as.double(xtx),
				xty        =as.double(xty),
				grid       =as.double(grid),
                workcpo    =as.double(workcpo),
				PACKAGE    ="DPpackage")


         #########################################################################################
         # save state
         #########################################################################################

           hpdf <- function(x)
           {
                alpha <- 0.05
                vec <- x
                n <- length(x)         
                alow <- rep(0,2)
                aupp <- rep(0,2)
                a <-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                           alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                return(c(a$alow[1],a$aupp[1]))
           }
    
           pdf <- function(x)
           {
                alpha <- 0.05
                vec <- x
                n <- length(x)         
                alow<-rep(0,2)
                aupp<-rep(0,2)
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                          alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
              return(c(a$alow[2],a$aupp[2]))
           }

           model.name <- "Bayesian semiparametric Rasch Model"		
                
           state <- list(alpha=foo$alpha,
	                     b=foo$b,
	                     bclus=foo$bclus,
	                     beta=foo$beta,
	                     mu=foo$mu,
	                     ncluster=foo$ncluster,
	                     sigma2=foo$sigma2,
                         ss=foo$ss)
         
           cpo <- matrix(foo$cpo,nrow=nsubject,ncol=p)
           cdfsave <- matrix(foo$cdfsave,nrow=nsave,ncol=ngrid)
           randsave <- matrix(foo$randsave,nrow=nsave,ncol=nsubject+1)
           thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=p+5)
 
           cdf.m <- apply(cdfsave,2,mean)
           cdf.l <- NULL
           cdf.u <- NULL
           if(compute.band)
           {
              limm <- apply(cdfsave, 2, hpdf)
              cdf.l <- limm[1,]
              cdf.u <- limm[2,]
           }

           pnames <- NULL
           for(i in 2:p)
           {
               pnames <- c(pnames,paste("beta",i,sep=""))
           }
           pnames <- c(pnames,"mean","variance","mu","sigma2","ncluster","alpha")
           colnames(thetasave) <- pnames
         
           qnames <- NULL
           for(i in 1:nsubject)
           {
               temp <- paste("theta(ID=",i,sep="")
               temp <- paste(temp,")",sep="")
               qnames <- c(qnames,temp)
           }
           qnames <- c(qnames,"theta(Prediction)")
           dimnames(randsave) <- list(NULL,qnames)
         
           coeff <- apply(thetasave, 2, mean)
         
           save.state <- list(thetasave=thetasave,randsave=randsave,cdfsave=cdfsave)

	       z <- list(call=cl,
                     y=y,
                     modelname=model.name,
                     cpo=cpo,
					 prior=prior,
                     mcmc=mcmc,
					 state=state,
                     save.state=save.state,
					 nrec=nrec,
                     nsubject=nsubject,
                     p=p,
                     alpharand=alpharand,
                     murand=murand,
                     sigmarand=sigmarand,
                     acrate=foo$acrate,
                     coefficients=coeff,
                     cdf=cdf.m,
					 cdf.l=cdf.l,
                     cdf.u=cdf.u,
					 grid=grid,
					 cpov=foo$cpov)
                 
           cat("\n\n")
 	       class(z)<-"DPrasch"
		   return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 10-08-2006.
###

"print.DPrasch"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) {
        cat("Posterior Inference of Parameters:\n")
        if(x$alpharand==1){
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)}
        if(x$alpharand==0){
        print.default(format(x$coefficients[1:(length(x$coefficients)-1)], digits = digits), print.gap = 2, 
            quote = FALSE)}

    }

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
   
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.DPrasch"<-function(object, hpd=TRUE, ...) 
{
    stde<-function(x)
    {
    	n<-length(x)
    	return(sd(x)/sqrt(n))
    }

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

    #nsave<-object$nsave
    #dimen<-length(object$coefficients)
    #thetasave<-matrix(object$save.state$thetasave,nrow=nsave, ncol=dimen)
    thetasave<-object$save.state$thetasave


### Difficulty parameters

    dimen1 <- object$p-1

    mat <- thetasave[,1:dimen1]

    coef.p <- object$coefficients[1:dimen1]
    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

    if(hpd)
    {             
         limm <- apply(mat, 2, hpdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }
    else
    {
         limm <- apply(mat, 2, pdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }

    names(coef.m) <- names(object$coefficients[1:dimen1])
    names(coef.sd) <- names(object$coefficients[1:dimen1])
    names(coef.se) <- names(object$coefficients[1:dimen1])
    names(coef.l) <- names(object$coefficients[1:dimen1])
    names(coef.u) <- names(object$coefficients[1:dimen1])

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd)
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


### Functionals

    mat <- thetasave[,(dimen1+1):(dimen1+2)]

    coef.p <- object$coefficients[(dimen1+1):(dimen1+2)]
    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

    if(hpd)
    {             
         limm <- apply(mat, 2, hpdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }
    else
    {
         limm <- apply(mat, 2, pdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }

    names(coef.m) <- names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.sd) <- names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.se) <- names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.l) <- names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.u) <- names(object$coefficients[(dimen1+1):(dimen1+2)])

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd)
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

    ans$functionals <- coef.table

    dimen1 <- dimen1 + 2 

### CPO
    ans$cpo<-object$cpo

### Baseline Information

    dimen2 <- object$murand+object$sigmarand
 
    if(dimen2>0)
    {
    
       if(dimen2==1)
       {
          if(object$murand==1)
          {
             mat<-matrix(thetasave[,dimen1+1],ncol=1) 
             coef.p<-object$coefficients[dimen1+1]
          }
          else
          {
             mat<-matrix(thetasave[,dimen1+2],ncol=1) 
             coef.p<-object$coefficients[dimen1+2]
          }
       }
       else
       {
          mat<-thetasave[,(dimen1+1):(dimen1+2)]
          coef.p<-object$coefficients[(dimen1+1):(dimen1+2)]
       }
    
       coef.m <-apply(mat, 2, median)    
       coef.sd<-apply(mat, 2, sd)
       coef.se<-apply(mat, 2, stde)

       if(hpd){             
            limm<-apply(mat, 2, hpdf)
            coef.l<-limm[1,]
            coef.u<-limm[2,]
       }
       else
       {
            limm<-apply(mat, 2, pdf)
            coef.l<-limm[1,]
            coef.u<-limm[2,]
       }

       coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
       if(hpd)
       {
            dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                   "95%HPD-Low","95%HPD-Upp"))
       }
       else
       {
            dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                   "95%CI-Low","95%CI-Upp"))
       }

       ans$base<-coef.table
    }

### Precision parameter

    if(is.null(object$prior$a0))
    {
      dimen3<-1
      coef.p<-object$coefficients[(dimen1+2+1)]
      mat<-matrix(thetasave[,(dimen1+2+1)],ncol=1)
    }
    else
    {
      dimen3<-2
      coef.p<-object$coefficients[(dimen1+2+1):(dimen1+2+2)]
      mat<-thetasave[,(dimen1+2+1):(dimen1+2+2)]

    }  

    coef.m <-apply(mat, 2, median)    
    coef.sd<-apply(mat, 2, sd)
    coef.se<-apply(mat, 2, stde)

    if(hpd){             
         limm<-apply(mat, 2, hpdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }
    else
    {
         limm<-apply(mat, 2, pdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }


    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    
    if(hpd)
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
    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject
    ans$acrate<-object$acrate

    class(ans) <- "summaryDPrasch"
    return(ans)
}


"print.summaryDPrasch"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nDifficulty parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$functionals)) {
        cat("\nFunctionals of the Random Effects:\n")
        print.default(format(x$functionals, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$base)) {
        cat("\nBaseline distribution:\n")
        print.default(format(x$base, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No baseline parameters\n")

    if (length(x$prec)) {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.DPrasch"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
{

fancydensplot1<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
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


   if(is(x, "DPrasch"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients[1:(x$p+1)]
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:n)
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$murand==1)
           {
               title1<-paste("Trace of","mu",sep=" ")
               title2<-paste("Density of","mu",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+2]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+2],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$sigmarand==1)
           {
               title1<-paste("Trace of","sigma2",sep=" ")
               title2<-paste("Density of","sigma2",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+3]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+3],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           title1<-paste("Trace of","ncluster",sep=" ")
           title2<-paste("Density of","ncluster",sep=" ")       
           plot(ts(x$save.state$thetasave[,x$p+4]),main=title1,xlab="MCMC scan",ylab=" ")
           hist(x$save.state$thetasave[,x$p+4],main=title2,xlab="values", ylab="probability",probability=TRUE)

           if(x$alpharand==1)
           {
               title1<-paste("Trace of","alpha",sep=" ")
               title2<-paste("Density of","alpha",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+5]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+5],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           title1 <- c("Density Estimate")
           title2 <- c("CDF Estimate")
           #fancydensplot1(x$save.state$randsave[,(x$nsubject+1)],hpd=hpd,main=title1,xlab="theta",col=col)
           plot(x$grid,x$cdf,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")
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
            if(poss==0 && param !="predictive")             
	        {
	           stop("This parameter is not present in the original model.\n")
	        }
	    
	        par(ask = ask)
	        layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))

	        if(param !="predictive")
	        {
                title1<-paste("Trace of",pnames[poss],sep=" ")
                title2<-paste("Density of",pnames[poss],sep=" ")       
                plot(ts(x$save.state$thetasave[,poss]),main=title1,xlab="MCMC scan",ylab=" ")
                if(param=="ncluster")
                {
                   hist(x$save.state$thetasave[,poss],main=title2,xlab="values", ylab="probability",probability=TRUE)
                }
                else
                {
                  fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
                }   
            }
            else
            {
               title1<-c("Density Estimate")
               title2<-c("CDF Estimate")
               #fancydensplot1(x$save.state$randsave[,(x$nsubject+1)],hpd=hpd,main=title1,xlab="theta",col=col)
               plot(x$grid,x$cdf,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")
            }                

        }
   }

}


