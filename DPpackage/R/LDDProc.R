### LDDProc.R                    
### Fit a linear dependent DP model for conditional ROC curve estimation.
###
### Copyright: Alejandro Jara, 2011-2012.
###
### Last modification: 12-02-2012.
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

"LDDProc"<-
function(y.d,z.d,y.nond,z.nond,zpred.d,zpred.nond=NULL,prior.d,prior.nond=NULL,mcmc,state,status,ngrid=100,grid=NULL,
         compute.band=FALSE,type.band="PD",data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("LDDProc")

"LDDProc.default"<-
function(y.d,
         z.d,
         y.nond,
         z.nond,
         zpred.d,
         zpred.nond=NULL,
         prior.d,
         prior.nond=NULL,
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
 
        state.nond <- state$state.nond
        state.d <- state$state.d

        if(is.null(zpred.nond))
        {
		   zpred.nond <-  zpred.d
        }

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
       # fitting the models
       #########################################################################################

         cat("\nFitting the model for the non-diseased population ","\n")

         if(is.null(prior.nond))
         {
            prior.nond <-  prior.d
         }


         prior.nond$rocc <- 1
         prior.nond$nroc <- ngrid

         fit.nond <- LDDPdensity(y.nond~z.nond-1,prior=prior.nond,mcmc=mcmc,state=state.nond,
                                 status=status,ngrid=ngrid,grid=grid,
                                 zpred=zpred.nond,compute.band=compute.band,
                                 type.band=type.band)


         cat("\nFitting the model for the diseased population ","\n")

         prior.d$rocc <- 2
         prior.d$nroc <- ngrid
         fit.d <- LDDPdensity(y.d~z.d-1,prior=prior.d,mcmc=mcmc,state=state.d,status=status,
							  ngrid=ngrid,grid=grid,zpred=zpred.d,compute.band=compute.band,
                              type.band=type.band)

       #########################################################################################
       # save state
       #########################################################################################

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

         if(!is.null(work.dir))
         {
            cat("\n\n Changing working directory back to ",old.dir,"\n")
            setwd(old.dir)
         }

         cpo.nond <- fit.nond$cpo
         cpo.d <- fit.d$cpo
         state.nond <- fit.nond$state
         state.d <- fit.d$state

         state <- list(state.nond=state.nond,state.d=state.d)

         model.name <-"Bayesian Conditional ROC Curve Estimation using LDDP Mixture of Normals"
         
         rocgrid <- fit.d$rocgrid
         rocp.m <- fit.d$rocp.m
         rocp.l <- fit.d$rocp.l
         rocp.h <- fit.d$rocp.h

         aucsave <- fit.d$save.state$aucsave

         aucp.m <- apply(aucsave,2,mean)
		 if(type.band!="HPD")
		 {
			limm <- apply(aucsave, 2, pdf)
			aucp.l <- limm[1,]
			aucp.h <- limm[2,]
		 }
         else
         {
			limm <- apply(aucsave, 2, hpdf)
			aucp.l <- limm[1,]
			aucp.h <- limm[2,]
		 }


		 z <- list(	modelname=model.name,
					call=cl,
					cpo.d=cpo.d,
                    cpo.nond=cpo.nond,
                    fit.d=fit.d,
                    fit.nond=fit.nond,
					prior.d=prior.d,
                    prior.nond=prior.nond,
					mcmc=mcmc,
                    rocgrid=rocgrid,
					aucp.m=aucp.m,
					aucp.l=aucp.l,
					aucp.h=aucp.h,
					rocp.m=rocp.m,
					rocp.l=rocp.l,
					rocp.h=rocp.h,
					state=state,
					work.dir=work.dir)

		 cat("\n\n")
		 class(z) <- c("LDDProc")
		 z 
}

###                    
### Tools
###
### Copyright: Alejandro Jara, 2011
### Last modification: 15-11-2011.
###


"print.LDDProc" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(c(log(x$cpo.d),log(x$cpo.nond))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    cat("\nPosterior Inference of Parameters - Non-diseased Group:\n")
    print.default(format(x$fit.nond$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nPosterior Inference of Parameters - Diseased Group:\n")
    print.default(format(x$fit.d$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nNumber of Observations in Non-diseased Group:",x$fit.nond$nrec)
    cat("\nNumber of Observations in Diseased Group:",x$fit.d$nrec)
    cat("\nNumber of Predictors:",x$fit.nond$p,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.LDDProc"<-function(object, hpd=TRUE, ...) 
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

    thetasave.nond <- object$fit.nond$save.state$thetasave
    thetasave.d <- object$fit.d$save.state$thetasave

    ans <- c(object[c("call", "modelname")])

    p <- object$fit.nond$p
    ans$p <- object$fit.nond$p

### CPO
    ans$cpo <-c(object$cpo.d,object$cpo.nond)

### Baseline Information - Nondiseased

    mat <- NULL
    coef.p <- NULL
    
    dimen1 <- p+p*(p+1)/2+1
    
    coef.p <- object$fit.nond$coefficients[1:dimen1]
    mat <- thetasave.nond[,1:dimen1]

    if(dimen1>0)
    {
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

       ans$base.nond <- coef.table
    }

### Precision parameter - Non-diseased

    dimen1 <- p+p*(p+1)/2+1
    
    if(is.null(object$prior$a0))
    {
      dimen2 <- 1
      coef.p <- object$fit.nond$coefficients[(dimen1+1)]
      mat <- matrix(thetasave.nond[,(dimen1+1)],ncol=1)
    }
    else
    {
      dimen2 <- 2
      coef.p <- object$fit.nond$coefficients[(dimen1+1):(dimen1+2)]
      mat <- thetasave.nond[,(dimen1+1):(dimen1+2)]

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

    ans$prec.nond <-coef.table


### Baseline Information - Diseased

    mat <- NULL
    coef.p <- NULL
    
    dimen1 <- p+p*(p+1)/2+1
    
    coef.p <- object$fit.d$coefficients[1:dimen1]
    mat <- thetasave.d[,1:dimen1]

    if(dimen1>0)
    {
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

       ans$base.d <- coef.table
    }

### Precision parameter - Diseased

    dimen1 <- p+p*(p+1)/2+1
    
    if(is.null(object$prior$a0))
    {
      dimen2 <- 1
      coef.p <- object$fit.d$coefficients[(dimen1+1)]
      mat <- matrix(thetasave.d[,(dimen1+1)],ncol=1)
    }
    else
    {
      dimen2 <- 2
      coef.p <- object$fit.d$coefficients[(dimen1+1):(dimen1+2)]
      mat <- thetasave.d[,(dimen1+1):(dimen1+2)]

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

    ans$prec.d <-coef.table


### Precision parameter - Diseased

    ans$nrec.nond <- object$fit.nond$nrec
    ans$nrec.d <- object$fit.d$nrec

    class(ans) <- "summaryLDDProc"
    return(ans)
}


"print.summaryLDDProc"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$base.nond)) {
        cat("\nBaseline distribution - Non-diseased Group:\n")
        print.default(format(x$base.nond, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No baseline parameters\n")


    if (length(x$base.d)) {
        cat("\nBaseline distribution - Diseased Group:\n")
        print.default(format(x$base.d, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No baseline parameters\n")

    if (length(x$prec.nond)) {
        cat("\nPrecision parameter - Non-diseased Group:\n")
        print.default(format(x$prec.nond, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$prec.d)) {
        cat("\nPrecision parameter - Diseased Group:\n")
        print.default(format(x$prec.d, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nNumber of Observations in Non-diseased Group:",x$nrec.nond)
    cat("\nNumber of Observations in Diseased Group:",x$nrec.d)
    cat("\nNumber of Predictors:",x$p,"\n")        
    cat("\n\n")
    invisible(x)
}


"plot.LDDProc" <- function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "LDDProc"))
   {
        if(is.null(param))
        {
           coef.p <- x$fit.nond$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:length(coef.p))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$fit.nond$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$fit.nond$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$fit.nond$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }

           coef.p <- x$fit.d$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:length(coef.p))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$fit.d$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$fit.d$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$fit.d$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }

           for(i in 1:x$fit.nond$npred)
           {
               if(x$fit.nond$compute.band)
               {
                  title1 <- paste("Density Prediction - Non-Diseased #",i,sep=" ")           
                  plot(x$fit.nond$grid,x$fit.nond$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
                  lines(x$fit.nond$grid,x$fit.nond$densp.l[i,],lty=2,lwd=2)
                  lines(x$fit.nond$grid,x$fit.nond$densp.m[i,],lty=1,lwd=3)

                  title1 <- paste("Density Prediction - Diseased #",i,sep=" ")           
                  plot(x$fit.d$grid,x$fit.d$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
                  lines(x$fit.d$grid,x$fit.d$densp.l[i,],lty=2,lwd=2)
                  lines(x$fit.d$grid,x$fit.d$densp.m[i,],lty=1,lwd=3)
			   }
               else
               {
				  title1 <- paste("Density Prediction - Non-Diseased #",i,sep=" ")           
				  plot(x$fit.nond$grid,x$fit.nond$densp.m[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")

				  title1 <- paste("Density Prediction - Diseased #",i,sep=" ")           
				  plot(x$fit.d$grid,x$fit.d$densp.m[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
			   }
           }


           for(i in 1:x$fit.nond$npred)
           {
               if(x$fit.nond$compute.band)
               {
                  title1 <- paste("ROC Curve Prediction #",i,sep=" ")           
                  plot(x$rocgrid,x$rocp.h[i,],main=title1,lty=2,type='l',lwd=2,
                       xlim=c(0,1),ylim=c(0,1),
                       xlab="False positive rate", 
                       ylab="True positive rate")
 
                  lines(x$rocgrid,x$rocp.l[i,],lty=2,lwd=2)
                  lines(x$rocgrid,x$rocp.m[i,],lty=1,lwd=3)

			   }
               else
               {
                  title1 <- paste("ROC Curve Prediction #",i,sep=" ")           
                  plot(x$rocgrid,x$rocp.m[i,],main=title1,lty=1,type='l',lwd=3,
                       xlim=c(0,1),ylim=c(0,1),
                       xlab="False positive rate", 
                       ylab="True positive rate")

			   }
           }

           
        }

        else
        {
            coef.p <- x$fit.nond$coefficients
			n<-length(coef.p)
			pnames<-names(coef.p)
			poss.nond<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss.nond=i
            }
            if(poss.nond==0 && param !="predictive")             
	        {
	           stop("This parameter is not present in the original model for Nondiseased.\n")
			}

            coef.p <- x$fit.d$coefficients
			n<-length(coef.p)
			pnames<-names(coef.p)
			poss.d<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss.d=i
            }
            if(poss.d==0 && param !="predictive")             
	        {
	           stop("This parameter is not present in the original model for Diseased.\n")
			}

			par(ask = ask)
	        layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))

	        if(param !="predictive")
	        {
                title1 <- paste("Trace of",pnames[poss.nond],sep=" ")
                title2 <- paste("Density of",pnames[poss.nond],sep=" ")       
                plot(ts(x$fit.nond$save.state$thetasave[,poss.nond]),main=title1,xlab="MCMC scan",ylab=" ")
                if(param=="ncluster")
                {
                   hist(x$fit.nond$save.state$thetasave[,poss.nond],main=title2,xlab="values", ylab="probability",probability=TRUE)
                }
                else
                {
                  fancydensplot1(x$fit.nond$save.state$thetasave[,poss.nond],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
                }   

                title1 <- paste("Trace of",pnames[poss.d],sep=" ")
                title2 <- paste("Density of",pnames[poss.d],sep=" ")       
                plot(ts(x$fit.nond$save.state$thetasave[,poss.d]),main=title1,xlab="MCMC scan",ylab=" ")
                if(param=="ncluster")
                {
                   hist(x$fit.nond$save.state$thetasave[,poss.d],main=title2,xlab="values", ylab="probability",probability=TRUE)
                }
                else
                {
                  fancydensplot1(x$fit.nond$save.state$thetasave[,poss.d],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
                }   

            }
            else
            {

				for(i in 1:x$fit.nond$npred)
				{
					if(x$fit.nond$compute.band)
					{
						title1 <- paste("Density Prediction - Non-Diseased #",i,sep=" ")           
						plot(x$fit.nond$grid,x$fit.nond$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
						lines(x$fit.nond$grid,x$fit.nond$densp.l[i,],lty=2,lwd=2)
						lines(x$fit.nond$grid,x$fit.nond$densp.m[i,],lty=1,lwd=3)

						title1 <- paste("Density Prediction - Diseased #",i,sep=" ")           
						plot(x$fit.d$grid,x$fit.d$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
						lines(x$fit.d$grid,x$fit.d$densp.l[i,],lty=2,lwd=2)
						lines(x$fit.d$grid,x$fit.d$densp.m[i,],lty=1,lwd=3)
					}
					else
					{
						title1 <- paste("Density Prediction - Non-Diseased #",i,sep=" ")           
						plot(x$fit.nond$grid,x$fit.nond$densp.m[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")

						title1 <- paste("Density Prediction - Diseased #",i,sep=" ")           
						plot(x$fit.d$grid,x$fit.d$densp.m[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
					}
				}


				for(i in 1:x$fit.nond$npred)
				{
					if(x$fit.nond$compute.band)
					{
						title1 <- paste("ROC Curve Prediction #",i,sep=" ")           
						plot(x$rocgrid,x$rocp.h[i,],main=title1,lty=2,type='l',lwd=2,
							xlim=c(0,1),ylim=c(0,1),
							xlab="False positive rate", 
							ylab="True positive rate")
 
						lines(x$rocgrid,x$rocp.l[i,],lty=2,lwd=2)
						lines(x$rocgrid,x$rocp.m[i,],lty=1,lwd=3)

					}
					else
					{
						title1 <- paste("ROC Curve Prediction #",i,sep=" ")           
						plot(x$rocgrid,x$rocp.m[i,],main=title1,lty=1,type='l',lwd=3,
							xlim=c(0,1),ylim=c(0,1),
							xlab="False positive rate", 
							ylab="True positive rate")

					}
				}

            }                
        }
   }

}



