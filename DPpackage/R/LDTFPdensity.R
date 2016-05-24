### LDTFPdensity.R                    
### Fit a linear dependent TF process for density regression.
###
### Copyright: Alejandro Jara, 2011-2012.
### Last modification: 11-11-2011.
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

"LDTFPdensity" <-
function(y,x,xtf,prediction,prior,mcmc,state,status,ngrid=100,grid=NULL,compute.band=FALSE,type.band="PD",
data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("LDTFPdensity")

"LDTFPdensity.default" <-
function(y,
         x,
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
         xce <- x
         pce <- ncol(xce)
         ptf <- ncol(xtf)

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

         xcepred <- prediction$xdenpred
         xtfpred <- prediction$xtfdenpred
         npredden <- nrow(xcepred)

         xcepredm <- prediction$xmedpred
         xtfpredm <- prediction$xtfmedpred
         npredmed <- nrow(xcepredm)

         quans <- prediction$quans
         if(is.null(quans)) quans <- c(0.03,0.50,0.97)

         if(is.null(grid))
         {
			miny <- min(y)
			maxy <- max(y)
			sdy <- sqrt(var(y))
			grid <- seq(miny-0.01*sdy,maxy+0.01*sdy,len=ngrid)
		 }
		 else
		 {
			grid <- as.vector(grid)
			ngrid <- length(grid)
		 }

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
       # MLE analysis
       #########################################################################################

          fit0 <- lm(formula = y ~ xce-1)
		  betace <- coefficients(fit0)
          cgkvar <- vcov(fit0)
          sigma2 <- summary(fit0)$sigma^2
         
       #########################################################################################
       # Prior information
       #########################################################################################

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

         betacepm <- prior$mub 
         precce <- solve(prior$Sb)

		 if(is.null(prior$tau1))
		 {
			tau <- c(-1,-1)
		 }
		 else
		 {
            tau <- c(prior$tau1,prior$tau2)
		 }

         tfprior <- prior$tfprior
         if(is.null(tfprior))tfprior <- 1

         if(tfprior==1)
         {  
	        gprior <- 2*nrec*solve(t(xtf)%*%xtf)
         }
         if(tfprior==2)
         {
			gprior <- 2*(1/nrec)*t(xtf)%*%xtf
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

         cpo <- matrix(0,nrow=nrec,ncol=2)
         densm <- matrix(0,nrow=npredden,ncol=ngrid)
         densl <- matrix(0,nrow=npredden,ncol=ngrid)
         densu <- matrix(0,nrow=npredden,ncol=ngrid)
         survmm <- matrix(0,nrow=npredden,ncol=ngrid)
         survml <- matrix(0,nrow=npredden,ncol=ngrid)
         survmu <- matrix(0,nrow=npredden,ncol=ngrid)

         qmm <- matrix(0,nrow=npredmed,3)
         qml <- matrix(0,nrow=npredmed,3)
		 qmu <- matrix(0,nrow=npredmed,3)

		 thetasave <- matrix(0, nrow=nsave, ncol=(pce+2))
         randsave <- matrix(0, nrow=nsave, ncol=((ntlr-1)*ptf))

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
	       betace <- state$betace
	       sigma2 <- state$sigma2
	       betatf <- state$betatf
		 }

       #########################################################################################
       # working space
       #########################################################################################

	     seed <- c(sample(1:29000,1),sample(1:29000,1))

         iflag <- rep(0,pce)
         iflagtf <- rep(0,ptf)
         nobsbc <- rep(0,ntprob)
         obsbc <- matrix(0,nrow=ntprob,ncol=nrec)
         c0 <- matrix(0,nrow=ptf,ncol=ptf)

         workm1 <- matrix(0,nrow=pce,ncol=pce)
         workvh1 <- rep(0,(pce*(pce+1)/2))
         workv1 <- rep(0,pce)

         worksam <- rep(0,nsave)
         worksam2 <- matrix(0,nrow=nsave,ncol=ngrid)
         worksam3 <- matrix(0,nrow=nsave,ncol=npredmed)

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

         foo <- .Fortran("ldtfpdensity",
                          nrec        = as.integer(nrec),
                          ptf         = as.integer(ptf),
                          pce         = as.integer(pce),
                          xtf         = as.double(xtf),
                          xce         = as.double(xce),
                          y           = as.double(y),
                          ngrid       = as.integer(ngrid),
                          npredden    = as.integer(npredden),
                          npredmed    = as.integer(npredmed),
                          grid        = as.double(grid),
                          xtfpred     = as.double(xtfpred),
                          xcepred     = as.double(xcepred),
						  xtfpredm    = as.double(xtfpredm),
                          xcepredm    = as.double(xcepredm),
                          quans       = as.double(quans),
                          maxm        = as.integer(maxm),
                          ntprob      = as.integer(ntprob),
                          ntlr        = as.integer(ntlr),
                          a0b0        = as.double(a0b0),
                          betacepm    = as.double(betacepm),
                          gprior      = as.double(gprior),
                          precce      = as.double(precce),
                          tau         = as.double(tau),
                          alpha       = as.double(alpha),
                          betace      = as.double(betace),
                          betatf      = as.double(betatf),
                          sigma2      = as.double(sigma2),
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
                          survmm      = as.double(survmm),
                          survml      = as.double(survml),
                          survmu      = as.double(survmu),
                          thetasave   = as.double(thetasave),
                          randsave    = as.double(randsave),
                          iflag       = as.integer(iflag),
                          iflagtf     = as.integer(iflagtf),
                          nobsbc      = as.integer(nobsbc), 
                          obsbc       = as.integer(obsbc),
                          c0          = as.double(c0),
                          workm1      = as.double(workm1),
                          workvh1     = as.double(workvh1),
                          workv1      = as.double(workv1),
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
  
       #########################################################################################
       # save state
       #########################################################################################

         if(!is.null(work.dir))
         {
            cat("\n\n Changing working directory back to ",old.dir,"\n")
            setwd(old.dir)
         }

         model.name <- "Linear dependent TF process model for density regression"
         
         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo <- cpom[,1]
         fso <- cpom[,2]
       
	     densm <- matrix(foo$densm,nrow=npredden,ncol=ngrid)
         densl <- NULL
         densu <- NULL

         qmm <- matrix(foo$qmm,nrow=npredmed,3)
         qml <- NULL
         qmu <- NULL

         survmm <- matrix(foo$survmm,nrow=npredden,ncol=ngrid)
         survml <- NULL
         survmu <- NULL

         if(compute.band)
         {
			densl <- matrix(foo$densl,nrow=npredden,ncol=ngrid)
			densu <- matrix(foo$densu,nrow=npredden,ncol=ngrid)

			qml <- matrix(foo$qml,nrow=npredmed,3)
			qmu <- matrix(foo$qmu,nrow=npredmed,3)

			survml <- matrix(foo$survml,nrow=npredden,ncol=ngrid)
			survmu <- matrix(foo$survmu,nrow=npredden,ncol=ngrid)
         }

 	     thetasave <- matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(pce+2))
         randsave <- matrix(foo$randsave,nrow=mcmc$nsave, ncol=((ntlr-1)*ptf))
 	 
         colnames(thetasave) <- c(colnames(xce),"sigma2","alpha")
         coeff <- apply(thetasave,2,mean)

         colnames(randsave) <- rep(colnames(xtf),(ntlr-1))

	     state <- list(alpha=foo$alpha,
	                   betace=foo$betace,
	                   sigma2=foo$sigma2,
	                   betatf=matrix(foo$betatf,nrow=ntlr,ncol=ptf),
                       nobsbc=foo$nobsbc,
                       obsbc=matrix(foo$obsbc,nrow=ntprob,ncol=nrec))

   	     save.state <- list(thetasave=thetasave,
                            randsave=randsave)


	     z <- list(modelname=model.name,
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
                   pce=foo$pce,
                   ptf=foo$ptf,
                   y=y,
                   x=xce,
                   xtf=xtf,
                   ngrid=ngrid,
                   npredden=npredden,
                   npredmed=npredmed,
                   grid=grid,
                   densm=densm,
                   densl=densl,
                   densu=densu,
                   qmm=qmm,
                   qml=qml,
                   qmu=qmu,
                   survmm=survmm,
                   survml=survml,
                   survmu=survmu)

	    cat("\n\n")
	    class(z) <- c("LDTFPdensity")
	    z 
}


###
### Tools for LDTFPdensity: print, summary, plot
###
### Copyright: Alejandro Jara, 2011
### Last modification: 11-11-2011.



"print.LDTFPdensity" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors for the Median:",x$pce)    
    cat("\nNumber of Predictors for the Tailfree Probabilities:",x$ptf,"\n")    
    cat("\n\n")
    invisible(x)
}


"plot.LDTFPdensity"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "LDTFPdensity"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:length(coef.p))
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
			   fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           for(i in 1:x$npredden)
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
			for(i in 1:x$npredden)
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




"summary.LDTFPdensity" <- function(object, hpd=TRUE, ...) 
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
 
    dimen1 <- object$pce
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


    ans$nrec <- object$nrec
    ans$pce <- object$pce
	ans$ptf <- object$ptf

    class(ans) <- "summaryLDTFPdensity"
    return(ans)
}


"print.summaryLDTFPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors for the Median:",x$pce)    
    cat("\nNumber of Predictors for the Tailfree Probabilities:",x$pce,"\n")    
    cat("\n\n")
    invisible(x)
}



