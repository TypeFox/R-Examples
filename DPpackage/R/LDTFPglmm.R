### LDTFPglmm.R                    
### Fit a linear dependent TF process for the random effects distribution
### in the context of a GLMM.
###
### Copyright: Alejandro Jara, 2011 - 2012.
### Last modification: 29-06-2012.
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
### The authors' contact information:
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

"LDTFPglmm" <-
function(y,x,roffset=NULL,g,family,xtf,prediction,prior,mcmc,state,status,ngrid=100,grid=NULL,compute.band=FALSE,type.band="PD",
data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("LDTFPglmm")

"LDTFPglmm.default" <-
function(y,
         x,
         roffset=NULL,
         g,
         family,
         xtf,
         prediction,
         prior,
         mcmc,
         state,
         status, 
		 ngrid=100,
		 grid=NULL,
		 compute.band=FALSE,
		 type.band="PD",
         data=sys.frame(sys.parent()),
         na.action=na.fail,
         work.dir=NULL)
{
       #########################################################################################
       # call parameters
       #########################################################################################

         m <- mcall <- cl <- match.call()

       #########################################################################################
       # data structure
       #########################################################################################

         nrec <- length(y)
         subject <- g
         nsubject <- length(table(subject))

         nrecx <- nrow(x)
         nrecg <- length(subject)

         p <- ncol(x)
         ptf <- ncol(xtf)

         nsubjectx <- nrow(xtf)

		 if(nrec != nrecx)
		 {
             stop("The response vector must be of the same length than the number of rows of x.\n")      
		 }

		 if(nrec != nrecg)
		 {
             stop("The response vector must be of the same length than the vector of group indicators.\n")      
		 }

		 if(nsubject != nsubjectx)
		 {
             stop("The number of groups must be equal than the number of rows of xtf.\n")      
		 }

         freqsub <- table(subject)
         maxni <- max(freqsub)
         idrec <- seq(1,nrec)
         datastr <- matrix(0,nrow=nsubject,ncol=maxni+1)
         datastr[,1] <- freqsub
         for(i in 1:nsubject)
         {
             for(j in 1:freqsub[i])
             {
                 datastr[i,(j+1)] <- idrec[subject==i][j] 
             }
         }

         if(is.null(roffset)) roffset <- rep(0,nrec)

       #########################################################################################
       # change working directory (if requested..)
       #########################################################################################
         if(!is.null(work.dir))
         {
            cat("\n Changing working directory to ",work.dir,"\n")
            old.dir <- getwd()  # by default work in current working directory
            setwd(work.dir)
         }

       #########################################################################################
       # prediction
       #########################################################################################

         xpred <- prediction$xpred
         xtfpred <- prediction$xtfpred

         npred <- nrow(xpred)
         px <- ncol(xpred)
         ptfx <- ncol(xtfpred)

		 if(p != px)
		 {
             stop("The number of columns in 'x' and 'xpred' must be the same.\n")      
		 }

		 if(ptf != ptfx)
		 {
             stop("The number of columns in 'xtf' and 'xtfpred' must be the same.\n")      
		 }

         quans <- prediction$quans
         if(is.null(quans)) quans <- c(0.03,0.50,0.97)

		 cband <- 0
		 if(compute.band)
		 {
			cband <- 1
		 }
			
		 tband <- 0
		 if(type.band=="HPD")
		 {
			tband <- 1
		 }

       #########################################################################################
       # PQL analysis
       #########################################################################################

		 library(nlme)
		 library(MASS)
         if(is.null(roffset))
         { 
            fit0 <- glmmPQL(fixed=y~x-1, random=~1 | subject, family=family, verbose = FALSE) 
         }
         else
         {
            fit0 <- glmmPQL(fixed=y~x-1 + offset(roffset), random=~1 | subject, family=family, verbose = FALSE) 
         }
         beta <- fit0$coefficients$fixed
         b <- as.vector(fit0$coefficients$random$subject)
         sigma2b <- getVarCov(fit0)[1,1]

         if(is.null(grid))
         {
            theta <- beta[1] + b
			miny <- min(theta)
			maxy <- max(theta)
			sdy <- sqrt(var(theta))
			grid <- seq(miny-0.01*sdy,maxy+0.01*sdy,len=ngrid)
		 }
		 else
		 {
			grid <- as.vector(grid)
			ngrid <- length(grid)
		 }

       #########################################################################################
       # Prior information
       #########################################################################################

       # Fixed effects

		 betapm <- prior$mub 
         prec <- solve(prior$Sb)

         sbt <- prec%*%betapm
         sb <- cbind(sbt,betapm)

         px <- nrow(prec)
 		 if(p != px)
		 {
             stop("The dimension of the prior covariance matrix for the fixed effects is not proper.\n")      
		 }

         px <- length(betapm)
 		 if(p != px)
		 {
             stop("The dimension of the prior mean vector for the fixed effects is not proper.\n")      
		 }

       # LTFP parameters

         maxm <- prior$maxm
         ntprob <- 0
         ntlr <- 0
         for(i in 1:maxm)
         {
             ntprob <- ntprob + 2**i
             ntlr <- ntlr +2**(i-1)
         }
         
	     if(is.null(prior$a0))
	     {
		    a0b0 <- c(-1,-1)
		    alpha <- prior$alpha
	     }
	     else
	     {
	 	    a0b0 <- c(prior$a0,prior$b0)
	 	    alpha <- 1
	     }

		 if(is.null(prior$taub1))
		 {
			taub <- c(-1,-1)
		 }
		 else
		 {
            taub <- c(prior$taub1,prior$taub2)
		 }

         tfprior <- prior$tfprior
         if(is.null(tfprior))tfprior <- 1

         if(tfprior==1)
         {  
	        gprior <- 2*nsubject*solve(t(xtf)%*%xtf)
         }
         if(tfprior==2)
         {
			gprior <- 2*(1/nsubject)*t(xtf)%*%xtf
         }
         if(tfprior==3)
         {
            gprior <- diag(1000,ptf)
		 }

       #########################################################################################
       # mcmc specification
       #########################################################################################

         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,cband,tband)
         nsave <- mcmc$nsave

       #########################################################################################
       # output
       #########################################################################################

         acrate <- 0 
         cpo <- matrix(0,nrow=nrec,ncol=2)
         densm <- matrix(0,nrow=npred,ncol=ngrid)
         densl <- matrix(0,nrow=npred,ncol=ngrid)
         densu <- matrix(0,nrow=npred,ncol=ngrid)

         qmm <- matrix(0,nrow=npred,3)
         qml <- matrix(0,nrow=npred,3)
		 qmu <- matrix(0,nrow=npred,3)

		 thetasave <- matrix(0, nrow=nsave, ncol=(p+2))
         randsave <- matrix(0, nrow=nsave, ncol=nsubject)
         tfpsave <- matrix(0, nrow=nsave, ncol=((ntlr-1)*ptf))

       #########################################################################################
       # parameters depending on status
       #########################################################################################

       	 if(status)
	     { 
			betatf <- matrix(0,nrow=ntlr,ncol=ptf)
	     }	
      	 else
	     {
	       alpha <- state$alpha
           b <- state$b
	       beta <- state$beta
	       sigma2b <- state$sigma2b
	       betatf <- state$betatf
		 }

       #########################################################################################
       # working space
       #########################################################################################

	     seed <- c(sample(1:29000,1),sample(1:29000,1))

         iflagtf <- rep(0,ptf)
         nobsbc <- rep(0,ntprob)
         obsbc <- matrix(0,nrow=ntprob,ncol=nsubject)
         c0 <- matrix(0,nrow=ptf,ncol=ptf)

         iflag <- rep(0,p)
         betac <- rep(0,p)
         workm1 <- matrix(0,nrow=p,ncol=p)
         workmp1 <- matrix(0,nrow=p,ncol=p)
         workmp2 <- matrix(0,nrow=p,ncol=p)
         workvh1 <- rep(0,(p*(p+1)/2))
         workv1 <- rep(0,p)
         workvp1 <- rep(0,p)

         worksam <- rep(0,nsave)
         worksam2 <- matrix(0,nrow=nsave,ncol=ngrid)
         worksam3 <- matrix(0,nrow=nsave,ncol=npred)
         fs <- rep(0,ngrid)

         workm2 <- matrix(0,nrow=ptf,ncol=ptf)
         workvh2 <- rep(0,(ptf*(ptf+1)/2))
         workv2 <- rep(0,ptf)
         workv3 <- rep(0,ptf)
         workv4 <- rep(0,ptf)
 
         k <- rep(0,maxm)
       
         prob <- rep(0,2**maxm)
         probc <- rep(0,2**maxm)


       #########################################################################################
       # calling the fortran code
       #########################################################################################

         if(family$family=="poisson")
         {
            if(family$link=="log")
            {

                as <- c(alpha,sigma2b,acrate)
                z <- cbind(x,roffset)

                foo <- .Fortran("ldtfpglmmpois",
                          maxni       = as.integer(maxni),
                          nrec        = as.integer(nrec),
                          nsubject    = as.integer(nsubject),
                          p           = as.integer(p),
                          subject     = as.integer(subject),
                          datastr     = as.integer(datastr),
                          x           = as.double(z),
                          y           = as.integer(y),
                          ptf         = as.integer(ptf),
                          xtf         = as.double(xtf),
                          ngrid       = as.integer(ngrid),
                          npred       = as.integer(npred),
                          grid        = as.double(grid),
                          xpred       = as.double(xpred),
                          xtfpred     = as.double(xtfpred),
                          quans       = as.double(quans),
                          prec        = as.double(prec),
                          sb          = as.double(sb),
                          maxm        = as.integer(maxm),
                          ntprob      = as.integer(ntprob),
                          ntlr        = as.integer(ntlr),
                          a0b0        = as.double(a0b0),
                          gprior      = as.double(gprior),
                          taub        = as.double(taub),
                          as          = as.double(as),
                          b           = as.double(b),
                          beta        = as.double(beta),
                          betatf      = as.double(betatf),
                          mcmc        = as.integer(mcmcvec),
                          nsave       = as.integer(nsave),
                          seed        = as.integer(seed),
                          cpo         = as.double(cpo),
                          densm       = as.double(densm),
                          densl       = as.double(densl),
                          densu       = as.double(densu),
                          qmm         = as.double(qmm), 
                          qml         = as.double(qml), 
                          qmu         = as.double(qmu), 
                          thetasave   = as.double(thetasave),
                          randsave    = as.double(randsave),
                          tfpsave     = as.double(tfpsave),
                          betac       = as.double(betac),
                          workmp1     = as.double(workmp1),
                          workmp2     = as.double(workmp2),
                          iflag       = as.integer(iflag),
                          iflagtf     = as.integer(iflagtf),
                          nobsbc      = as.integer(nobsbc), 
                          obsbc       = as.integer(obsbc),
                          c0          = as.double(c0),
                          workm1      = as.double(workm1),
                          workvh1     = as.double(workvh1),
                          workv1      = as.double(workv1),
                          workvp1     = as.double(workvp1),
                          worksam     = as.double(worksam),
                          worksam2    = as.double(worksam2),
						  worksam3    = as.double(worksam3),
                          fs          = as.double(fs),
                          workm2      = as.double(workm2),
                          workvh2     = as.double(workvh2),
                          workv2      = as.double(workv2),
                          workv3      = as.double(workv3),
                          workv4      = as.double(workv4),
                          k           = as.integer(k),
                          prob        = as.double(prob),
                          probc       = as.double(probc),
                          PACKAGE     = "DPpackage")	

                foo$alpha <- foo$as[1]
                foo$sigma2b <- foo$as[2]
                foo$acrate <- foo$as[3]
                acrate <- foo$acrate

            }
         }
  
       #########################################################################################
       # save state
       #########################################################################################

         if(!is.null(work.dir))
         {
            cat("\n\n Changing working directory back to ",old.dir,"\n")
            setwd(old.dir)
         }

         model.name <- "GLMM using a LDTFP prior for the random effects"

         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo <- cpom[,1]
         fso <- cpom[,2]
       
	     densm <- matrix(foo$densm,nrow=npred,ncol=ngrid)
         densl <- NULL
         densu <- NULL

         qmm <- matrix(foo$qmm,nrow=npred,3)
         qml <- NULL
         qmu <- NULL


         if(compute.band)
         {
			densl <- matrix(foo$densl,nrow=npred,ncol=ngrid)
			densu <- matrix(foo$densu,nrow=npred,ncol=ngrid)

			qml <- matrix(foo$qml,nrow=npred,3)
			qmu <- matrix(foo$qmu,nrow=npred,3)

         }

 	     thetasave <- matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+2))
         randsave <- matrix(foo$randsave,nrow=mcmc$nsave, ncol=nsubject)
         tfpsave <- matrix(foo$tfpsave,nrow=mcmc$nsave, ncol=((ntlr-1)*ptf))
 	 
         colnames(thetasave) <- c(colnames(x),"sigma2b","alpha")
         coeff <- apply(thetasave,2,mean)

         colnames(tfpsave) <- rep(colnames(xtf),(ntlr-1))

	     state <- list(alpha=foo$alpha,
	                   beta=foo$beta,
                       b=foo$b,
	                   sigma2b=foo$sigma2b,
	                   betatf=matrix(foo$betatf,nrow=ntlr,ncol=ptf),
                       nobsbc=foo$nobsbc,
                       obsbc=matrix(foo$obsbc,nrow=ntprob,ncol=nsubject))

   	     save.state <- list(thetasave=thetasave,
                            randsave=randsave,
                            tfpsave=tfpsave)


	     z <- list(modelname=model.name,
                   acrate=acrate,
				   coefficients=coeff,
	               call=cl,
                   compute.band=compute.band,
	               cpo=cpo,
	               fso=fso,
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   save.state=save.state,
                   nrec=foo$nrec,
                   nsubject=foo$nsubject,
                   p=foo$p,
                   ptf=foo$ptf,
                   y=y,
                   x=x,
                   xtf=xtf,
                   ngrid=ngrid,
                   npred=npred,
                   grid=grid,
                   densm=densm,
                   densl=densl,
                   densu=densu,
                   qmm=qmm,
                   qml=qml,
                   qmu=qmu)

	    cat("\n\n")
	    class(z) <- c("LDTFPglmm")
	    z 
}


###
### Tools for LDTFPglmm: print, summary, plot
###
### Copyright: Alejandro Jara, 2011
### Last modification: 11-11-2011.


"print.LDTFPglmm" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    

    cat("\nNumber of groups:",x$nsubject)
    cat("\nNumber of observations:",x$nrec)
    cat("\nNumber of predictors for the fixed effects:",x$p)    
    cat("\nNumber of predictors for the tailfree probabilities:",x$ptf,"\n")    
    cat("\n\n")
    invisible(x)
}


"plot.LDTFPglmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
{

fancydensplot1<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
# Author: AJV, 2007
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


   if(is(x, "LDTFPglmm"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           pp <- length(coef.p)
           if(is.null(x$prior$a0))
           {
              pp <- pp - 1
           }

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:pp)
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
			   fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           for(i in 1:x$npred)
           {
               if(x$compute.band)
               {
                  title1 <- paste("Density Prediction #",i,sep=" ")           
                  plot(x$grid,x$densu[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
                  lines(x$grid,x$densl[i,],lty=2,lwd=2)
                  lines(x$grid,x$densm[i,],lty=1,lwd=3)
			   }
               else
               {
				  title1 <- paste("Density Prediction #",i,sep=" ")           
				  plot(x$grid,x$densm[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
			   }
           }
           
        }

        else
        {
            coef.p <- x$coefficients
			n <- length(coef.p)
			pnames <- names(coef.p)
			poss <- 0 
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
                title1 <- paste("Trace of",pnames[poss],sep=" ")
                title2 <- paste("Density of",pnames[poss],sep=" ")       
                plot(ts(x$save.state$thetasave[,poss]),main=title1,xlab="MCMC scan",ylab=" ")
				fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
		}
		else
		{
			for(i in 1:x$npred)
			{
				if(x$compute.band)
				{
						title1 <- paste("Density Prediction #",i,sep=" ")           
						plot(x$grid,x$densu[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
						lines(x$grid,x$densl[i,],lty=2,lwd=2)
						lines(x$grid,x$densm[i,],lty=1,lwd=3)
				}
				else
				{
						title1 <- paste("Density Prediction #",i,sep=" ")           
						plot(x$grid,x$densm[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
				}
			}
		}                
	}
}

}


"summary.LDTFPglmm" <- function(object, hpd=TRUE, ...) 
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

    thetasave <- object$save.state$thetasave

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo <- object$cpo

### Median information
 
    dimen1 <- object$p
    if(dimen1==1)
    {
       mat <- matrix(thetasave[,1],ncol=1)
	}
    else
    {
       mat <- thetasave[,1:dimen1]
    }

    coef.p <- object$coefficients[1:dimen1]    
	coef.m <- apply(mat, 2, median)    
	coef.sd <- apply(mat, 2, sd)
	coef.se <- apply(mat, 2, stde)

	if(hpd){             
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

	ans$coeff <- coef.table


### Baseline Information

    mat <- matrix(thetasave[,(dimen1+1)],ncol=1)

    coef.p <- object$coefficients[(dimen1+1)]    
	coef.m <- apply(mat, 2, median)    
	coef.sd <- apply(mat, 2, sd)
	coef.se <- apply(mat, 2, stde)

	if(hpd){             
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

	ans$base <- coef.table

### Precision parameter

    if(is.null(object$prior$a0))
    {
       ans$prec <- NULL
	}
    else
	{

		mat <- matrix(thetasave[,(dimen1+2)],ncol=1)

        coef.p <- object$coefficients[(dimen1+2)]    
		coef.m <- apply(mat, 2, median)    
		coef.sd <- apply(mat, 2, sd)
		coef.se <- apply(mat, 2, stde)

		if(hpd){             
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

		ans$prec <- coef.table
	}

    ans$acrate <- object$acrate
    ans$nrec <- object$nrec
    ans$p <- object$p
	ans$ptf <- object$ptf

    class(ans) <- "summaryLDTFPglmm"
    return(ans)
}




"print.summaryLDTFPglmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    cat("\nPosterior Inference of Median Regression Parameters:\n")
    print.default(format(x$coeff, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nPosterior Inference of Baseline Variance:\n")
    print.default(format(x$base, digits = digits), print.gap = 2, 
            quote = FALSE)

    if (length(x$prec)) {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    

    cat("\nNumber of groups:",x$nsubject)
    cat("\nNumber of observations:",x$nrec)
    cat("\nNumber of predictors for the fixed effects:",x$p)    
    cat("\nNumber of predictors for the tailfree probabilities:",x$ptf,"\n")    
    cat("\n\n")
    invisible(x)
}

