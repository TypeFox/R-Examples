### LDDPsurvival.R                    
### Fit a linear dependent DP model for survival modelling.
###
### Copyright: Alejandro Jara, Peter Mueller and Gary Rosner, 2009-2012.
###
### Last modification: 15-07-2012.
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
###      Peter Mueller
###      Department of Mathematics
###      The University of Texas at Austin
###      1, University Station, C1200
###      Austin, TX 78712, USA
###      Voice: (512) 471-7168  URL  : http://www.ma.utexas.edu/users/pmueller/
###      Fax  : (512) 471-9038  Email: pmueller@math.utexas.edu
###
###      Gary L. Rosner
###      Division of Oncology Biostatistics/Bioinformatics
###      The Sidney Kimmel Comprehensive Cancer Center
###      Johns Hopkins
###      550 North Broadway, Suite 1103
###      Baltimore, Maryland  21205-2013
###      Voice: (410) 955-4884  URL  : http://www.hopkinskimmelcancercenter.org/index.cfm/cID/1686/mpage/expertdata.cfm/expID/593
###      Fax  : (410) 955-0859  Email: grosner@jhmi.edu
###


"LDDPsurvival" <-
function(formula,zpred,prior,mcmc,state,status,grid,data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("LDDPsurvival")

"LDDPsurvival.default" <-
function(formula,
         zpred,
         prior,
         mcmc,
         state,
         status, 
         grid,
         data=sys.frame(sys.parent()),
         na.action=na.fail,
         work.dir=NULL)
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
       # data structure
       #########################################################################################
		 ymat <- model.response(mf,"numeric")
         nrec <- nrow(ymat)
		 z <- model.matrix(formula)
         p <- ncol(z)

         typed <- rep(2,nrec)
         for(i in 1:nrec)
         {
             typed[ymat[,1]==-999] <- 1
             typed[ymat[,2]==-999] <- 3
             typed[ymat[,1]==ymat[,2]] <- 4
         }

       #########################################################################################
       # MLE analysis
       #########################################################################################

          ymat2 <- ymat
          ymat2[ymat2[,1]==-999,1] <- 0
          ymat2[ymat2[,2]==-999,2] <- 10^10
          ll <- ymat2[,1]
          rr <- ymat2[,2]
          library(survival)
          fit0 <- survreg(formula = Surv(time=ll,time2=rr,type="interval2") ~ z-1, dist = "lognormal")

		  betaw <- coefficients(fit0)
          cgkvar <- vcov(fit0)[1:p,1:p]
          sigma2w <- fit0$scale

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
         npred <- nrow(zpred)
         ngrid <- length(grid)

       #########################################################################################
       # Prior information
       #########################################################################################

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

         sbeta0i <- solve(prior$S0)
         m0 <- prior$m0

         tau1 <- prior$tau1
         taus1 <- prior$taus1
         taus2 <- prior$taus2
         
         nu <- prior$nu
         psiinv <- prior$psiinv 
 
       #########################################################################################
       # mcmc specification
       #########################################################################################

         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave <- mcmc$nsave

       #########################################################################################
       # Starting values
       #########################################################################################
         ncluster <- 1
         ss <- rep(1,nrec)

         y <- rep(0,nrec)
         y[typed==2] <- log((ymat[typed==2,1]+ymat[typed==2,2])/2)
         y[typed==1] <- log(ymat[typed==1,2]/2)
		 y[typed==3] <- log(ymat[typed==3,1]+1)
         y[typed==4] <- log(ymat[typed==4,1])

         betaclus <- matrix(0,nrow=nrec+100,ncol=p)
         sigmaclus <- rep(0,nrec+100)
         
         betaclus[1,] <- betaw
         sigmaclus[1] <- sigma2w

         sb <- cgkvar
         mub <- betaw
         tau2 <- 2*sigma2w

		 censor <- sum(typed!=4)		
		 if(censor>0)censor <- 1
         
       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2)

         denspm <- matrix(0,nrow=npred,ncol=ngrid)
         denspl <- matrix(0,nrow=npred,ncol=ngrid)
         densph <- matrix(0,nrow=npred,ncol=ngrid)

         hazpm <- matrix(0,nrow=npred,ncol=ngrid)
         hazpl <- matrix(0,nrow=npred,ncol=ngrid)
         hazph <- matrix(0,nrow=npred,ncol=ngrid)

         meanfpm <- rep(0,npred)
         meanfpl <- rep(0,npred)
         meanfph <- rep(0,npred)

         survpm <- matrix(0,nrow=npred,ncol=ngrid)
         survpl <- matrix(0,nrow=npred,ncol=ngrid)
         survph <- matrix(0,nrow=npred,ncol=ngrid)

         thetasave <- matrix(0,nrow=nsave,ncol=(p+(p*(p+1)/2)+3))
         randsave <- matrix(0,nrow=nsave,ncol=nrec*(p+1))
         survsave <- matrix(0,nrow=nsave,ncol=npred*ngrid)

       #########################################################################################
       # parameters depending on status
       #########################################################################################


         if(status==FALSE) 
         {
            alpha <- state$alpha
            ncluster <- state$ncluster
            ss <- state$ss
            betaclus <- state$betaclus
            sigmaclus <- state$sigmaclus
            tau2 <- state$tau2
            mub <- state$mub
            sb <- state$sb
            y <- state$y
         }    

       #########################################################################################
       # working space
       #########################################################################################
         
         cstrt <- matrix(0,nrow=nrec,ncol=nrec)
         ccluster <- rep(0,nrec)
         iflagp <- rep(0,p)
         betam <- rep(0,p)
         betawork <- rep(0,p)
         prob <- rep(0,(nrec+100))
         workmh1 <- rep(0,p*(p+1)/2)
         workmh2 <- rep(0,p*(p+1)/2)
         workv1 <- rep(0,p)
         xtx <- matrix(0,nrow=p,ncol=p)
         xtx2 <- matrix(0,nrow=p,ncol=p)
         xty <- rep(0,p)
         xty2 <- rep(0,p)
         seed <- c(sample(1:29000,1),sample(1:29000,1))

         fs <- rep(0,ngrid) 
         fm <- rep(0,npred)

         worksam <- rep(0,nsave)
         
         workcpo <- rep(0,nrec)
         workcpo2 <- rep(0,nrec)

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         foo <- .Fortran("lddpsurvival",
                          nrec      = as.integer(nrec),
                          p         = as.integer(p),
                          censor    = as.integer(censor),
                          typed     = as.integer(typed),
                          ymat      = as.double(ymat),
                          z         = as.double(z),
                          ngrid     = as.integer(ngrid),
                          npred     = as.integer(npred),
                          grid      = as.double(grid),
                          zpred     = as.double(zpred),
                          a0b0      = as.double(a0b0), 
                          tau1      = as.double(tau1),
                          taus1     = as.double(taus1),
                          taus2     = as.double(taus2),
                          m0        = as.double(m0),
                          sbeta0i   = as.double(sbeta0i),
                          nu        = as.integer(nu),
                          psiinv    = as.double(psiinv),
                          ncluster  = as.integer(ncluster), 
                          ss        = as.integer(ss),
                          alpha     = as.double(alpha),
                          betaclus  = as.double(betaclus),
                          sigmaclus = as.double(sigmaclus),
                          mub       = as.double(mub),
                          sb        = as.double(sb), 
                          tau2      = as.double(tau2),
                          y         = as.double(y),
                          cpo       = as.double(cpo),
                          thetasave = as.double(thetasave),
                          randsave  = as.double(randsave),
                          survsave  = as.double(survsave),
                          denspm    = as.double(denspm),
                          denspl    = as.double(denspl),
                          densph    = as.double(densph),
                          hazpm     = as.double(hazpm),
                          hazpl     = as.double(hazpl),
                          hazph     = as.double(hazph),
                          meanfpm   = as.double(meanfpm),
                          meanfpl   = as.double(meanfpl),
                          meanfph   = as.double(meanfph),
                          survpm    = as.double(survpm),
                          survpl    = as.double(survpl),
                          survph    = as.double(survph),
                          mcmc      = as.integer(mcmcvec),
                          nsave     = as.integer(nsave),
                          seed      = as.integer(seed),
                          cstrt     = as.integer(cstrt),
                          ccluster  = as.integer(ccluster),
                          iflagp    = as.integer(iflagp),
                          betam     = as.double(betam),
                          betawork  = as.double(betawork),
                          prob      = as.double(prob),
                          workmh1   = as.double(workmh1),
                          workmh2   = as.double(workmh2),
                          workv1    = as.double(workv1),
                          xtx       = as.double(xtx),
                          xtx2      = as.double(xtx2),
                          xty       = as.double(xty),
                          xty2      = as.double(xty2),
                          fs        = as.double(fs),
                          fm        = as.double(fm),
                          worksam   = as.double(worksam),
                          workcpo   = as.double(workcpo),
                          workcpo2  = as.double(workcpo2),
                          PACKAGE="DPpackage")	
         
       #########################################################################################
       # save state
       #########################################################################################

         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo <- cpom[,1]         
         fso <- cpom[,2]

         model.name <- "Bayesian semiparametric survival model using a LDDP mixture"
         
         state <- list(alpha=foo$alpha,
	               betaclus=matrix(foo$betaclus,nrow=nrec+100,ncol=p),
	               sigmaclus=foo$sigmaclus,
	               ss=foo$ss,
	               ncluster=foo$ncluster,
	               mub=foo$mub,
	               sb=matrix(foo$sb,nrow=p,ncol=p),
	               tau2=foo$tau2,
                   y=foo$y)

         denspm <- matrix(foo$denspm,nrow=npred,ncol=ngrid)
         denspl <- matrix(foo$denspl,nrow=npred,ncol=ngrid)
         densph <- matrix(foo$densph,nrow=npred,ncol=ngrid)

         hazpm <- matrix(foo$hazpm,nrow=npred,ncol=ngrid)
         hazpl <- matrix(foo$hazpl,nrow=npred,ncol=ngrid)
         hazph <- matrix(foo$hazph,nrow=npred,ncol=ngrid)

         meanfpm <- foo$meanfpm
         meanfpl <- foo$meanfpl
         meanfph <- foo$meanfph

         survpm <- matrix(foo$survpm,nrow=npred,ncol=ngrid)
         survpl <- matrix(foo$survpl,nrow=npred,ncol=ngrid)
         survph <- matrix(foo$survph,nrow=npred,ncol=ngrid)

         randsave <- matrix(foo$randsave,nrow=nsave,ncol=nrec*(p+1))
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=(p+(p*(p+1)/2)+3))
         survsave <- matrix(foo$survsave,nrow=nsave,ncol=(npred*ngrid))

         survpmed <- matrix(apply(survsave,2,median),nrow=npred,ncol=ngrid,byrow=TRUE)

         coeffname <- dimnames(z)[[2]]

         pnames1 <- NULL
         for(i in 1:p)
         {
             pnames1 <- c(pnames1,paste("mub",coeffname[i],sep=""))
         }
         
         pnames2 <- NULL
         for(i in 1:p)
         {
            for(j in i:p)
            {
                tmp <- paste("sb",coeffname[i],sep="")
                tmp <- paste(tmp,coeffname[j],sep=":")
                pnames2 <- c(pnames2,tmp)
            }    
         }

         pnames <- c(pnames1,pnames2,"tau2","ncluster","alpha")
         colnames(thetasave) <- pnames


         coeff <- apply(thetasave, 2, mean)

         renames <- NULL
         for(i in 1:nrec)
         {
             tmp <- paste(coeffname,i,sep=":")
             renames <- c(renames,tmp)
             tmp <- paste("sigma2",i,sep=":")
             renames <- c(renames,tmp)
         }
         colnames(randsave) <- renames
         
         save.state <- list(thetasave=thetasave,
                            randsave=randsave)

         z<-list(modelname=model.name,
	             call=cl,
	             cpo=cpo,
	             coefficients=coeff,
	             fso=fso,
                 prior=prior,
                 mcmc=mcmc,
                 nrec=foo$nrec,
                 p=foo$p,
                 z=z,
                 ngrid=ngrid,
		         npred=npred,
		         zpred=zpred,
		         grid=grid,
		         densp.m=denspm,
                 densp.l=denspl,
                 densp.h=densph,
				 hazp.m=hazpm,
                 hazp.l=hazpl,
                 hazp.h=hazph,
                 meanfp.m=meanfpm,
                 meanfp.l=meanfpl,
                 meanfp.h=meanfph,
		         survp.m=survpm,
                 survp.med=survpmed,
                 survp.l=survpl,
                 survp.h=survph,
                 state=state,
                 save.state=save.state)

	 cat("\n\n")
	 class(z) <- c("LDDPsurvival")
	 z 
}



###                    
### Tools
###
### Copyright: Alejandro Jara, 2008
### Last modification: 02-06-2008.
###


"print.LDDPsurvival" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Predictors:",x$p,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.LDDPsurvival"<-function(object, hpd=TRUE, ...) 
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

    thetasave<-object$save.state$thetasave

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo<-object$cpo

### Baseline Information

    mat<-NULL
    coef.p<-NULL
    
    dimen1 <- object$p+object$p*(object$p+1)/2+1
    
    coef.p <- object$coefficients[1:dimen1]
    mat <- thetasave[,1:dimen1]

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

       ans$base<-coef.table
    }

### Precision parameter

    dimen1<-object$p+object$p*(object$p+1)/2
    
    if(is.null(object$prior$a0))
    {
      dimen2<-1
      coef.p<-object$coefficients[(dimen1+1)]
      mat<-matrix(thetasave[,(dimen1+1)],ncol=1)
    }
    else
    {
      dimen2<-2
      coef.p<-object$coefficients[(dimen1+1):(dimen1+2)]
      mat<-thetasave[,(dimen1+1):(dimen1+2)]

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
    ans$p<-object$p

    class(ans) <- "summaryLDDPsurvival"
    return(ans)
}


"print.summaryLDDPsurvival"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors:",x$p,"\n")        
    cat("\n\n")
    invisible(x)
}



"plot.LDDPsurvival"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "LDDPsurvival"))
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
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }

           for(i in 1:x$npred)
           {
               title1 <- paste("Density Prediction #",i,sep=" ")           
               plot(x$grid,x$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
               lines(x$grid,x$densp.l[i,],lty=2,lwd=2)
               lines(x$grid,x$densp.m[i,],lty=1,lwd=3)
           }

           for(i in 1:x$npred)
           {
               title1 <- paste("Survival Prediction #",i,sep=" ")           
               plot(x$grid,x$survp.h[i,],main=title1,type="l",xlab="time",ylab="survival",lty=2,lwd=2,ylim=c(0,1))
               lines(x$grid,x$survp.l[i,],lty=2,lwd=2)
               lines(x$grid,x$survp.m[i,],lty=1,lwd=3)
           }

           for(i in 1:x$npred)
           {
               title1 <- paste("Hazard Prediction #",i,sep=" ")           
               plot(x$grid,x$hazp.h[i,],main=title1,type="l",xlab="time",ylab="hazard",lty=2,lwd=2)
               lines(x$grid,x$hazp.l[i,],lty=2,lwd=2)
               lines(x$grid,x$hazp.m[i,],lty=1,lwd=3)
           }
           
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
               for(i in 1:x$npred)
               {
                   title1 <- paste("Density Prediction #",i,sep=" ")           
                   plot(x$grid,x$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
                   lines(x$grid,x$densp.l[i,],lty=2,lwd=2)
                   lines(x$grid,x$densp.m[i,],lty=1,lwd=3)
               }
            }                
        }
   }

}




