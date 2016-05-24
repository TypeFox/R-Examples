### DPcdencity.R                   
### Fit a Dirichlet Process Mixture of Normals model for conditional density
### estimation.
###
### Copyright: Alejandro Jara, 2008-2012.
###
### Last modification: 27-08-2010.
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


"DPcdensity"<-
function(y,x,xpred,ngrid=100,grid=NULL,compute.band=FALSE,type.band="PD",prior,mcmc,state,status,data=sys.frame(sys.parent()),work.dir=NULL)
UseMethod("DPcdensity")

DPcdensity.default<-
function(y,
         x,
         xpred,
         ngrid=100,
		 grid=NULL,
		 compute.band=FALSE,
		 type.band="PD",
         prior,
         mcmc,
         state,
         status,
         data=sys.frame(sys.parent()),
         work.dir=NULL)
{
       #########################################################################################
       # call parameters
       #########################################################################################
         cl <- match.call()

       #########################################################################################
       # response
       #########################################################################################
         nrec <- length(y)
          
       #########################################################################################
       # add covariates to the data matrix
       #########################################################################################
         x <- as.matrix(x)
         nx <- dim(x)[2]
         z <- cbind(y,x)
         nvar <- nx + 1  

       #########################################################################################
       # identify missing values to be imputed
       #########################################################################################
         nmiss <- sum(is.na(z))
         nmissi <- 1
         missp <- NULL
         for(i in 1:nrec)
         {
             nmi <- sum(is.na(z[i,]))
             nmi2 <- seq(1,nvar)[is.na(z[i,])] 
             if(nmi>0)
             {
                for(j in 1:nmi)
                { 
                    missp <- rbind(missp,c(i,nmi2[j]))  
                }
             }  
         }
         
         if(nmiss==0)
         {
            nmissi <- 0
            nmiss <- 1
            missp <- matrix(0,nrow=nmiss,2)
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
       # prediction
       #########################################################################################
         xpred <- as.matrix(xpred)
         npred <- nrow(xpred)
         pz <- ncol(xpred)  

         if(pz != nx)
         {
            cat("\n *** Error: ncol(xpred) != ncol(x). Need to match.\n")
            return(-1)
         }
 
		 if(is.null(grid))
		 { 
			yy <- na.omit(y)  
			miny <- min(yy)
			maxy <- max(yy)
			vary <- var(yy)
			grid <- seq(from=miny-0.25*sqrt(vary),to=maxy+0.25*sqrt(vary),length.out=ngrid)
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
			
		 tband <- 1
		 if(type.band!="HPD")
		 {
			tband <- 2
		 }


       #########################################################################################
       # prior information
       #########################################################################################

         if(is.null(prior$a0))      
		 {
            a0 <--1
			b0 <--1 
            alpha<-prior$alpha
            alpharand<-0
         }
         else
         {
            a0 <- prior$a0
            b0 <- prior$b0
			alpha <- 1
			alpharand <- 1
		 }
         a0b0 <- c(a0,b0)
  	 
		 if(is.null(prior$nu2))
		 {
            psiinv1 <- matrix(prior$psiinv1,nvar,nvar)
            psiinv2 <- psiinv1
            psi1 <- matrix(solve(psiinv1),nvar,nvar)
            nuvec <- c(prior$nu1,-1)
            psi1rand <- 0
		 }
		 else
		 {
			psiinv1 <- matrix(var(z),nvar,nvar)
            psi1 <- matrix(solve(psiinv1),nvar,nvar)
            psiinv2 <- matrix(prior$psiinv2,nvar,nvar)
			nuvec <- c(prior$nu1,prior$nu2)
			psi1rand <- 1
		 }
  	 
		if(is.null(prior$m2) && is.null(prior$s2))
		{
			s2inv <- matrix(0,nrow=nvar,ncol=nvar) 
			s2invm2 <- matrix(0,nrow=nvar,ncol=1)
			m1 <- prior$m1
			m1rand <- 0
		}
		else
		{
            s2inv <- solve(prior$s2)
            s2invm2 <- s2inv%*%prior$m2
			m1 <- rep(0,nvar)
            for(i in 1:nvar)
      	    {
                m1[i] <- mean(z[,i])+rnorm(1,0,100)
       	    }
            m1rand <- 1
         }     

         if(is.null(prior$tau1) && is.null(prior$tau2)) 
         {
            tau <- c(-2,-2)
            k0 <- prior$k0
            k0rand <- 0
         }
         else
         {
            tau <- c(prior$tau1,prior$tau2)
            k0 <- rgamma(1,shape=prior$tau1,scale=prior$tau2)
            k0rand <- 1
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
         thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
         denspm <- matrix(0,nrow=npred,ncol=ngrid)
         denspl <- matrix(0,nrow=npred,ncol=ngrid)
         densph <- matrix(0,nrow=npred,ncol=ngrid)

         meanfpm <-rep(0,npred) 
         meanfpl <-rep(0,npred) 
         meanfph <-rep(0,npred) 

       #########################################################################################
       # parameters depending on status
       #########################################################################################
         nuniq <- nvar*(nvar+1)/2

         if(status==TRUE)
		 {
            muclus <- matrix(0,nrow=nrec+100,ncol=nvar)
            sigmaclus <- matrix(0,nrow=nrec+100,ncol=nuniq)
            for(i in 1:1)
            {
                counter <- 0
                for(j in 1:nvar)
                {
                    muclus[i,j] <- m1[j]
                    for(k in j:nvar)
                    {
          		counter<-counter+1
          		sigmaclus[i,counter]<-psiinv1[j,k]
                    }
                }
             }
             ncluster <- 1
             ss <- rep(1,nrec)
		 }
	 
      	 if(status==FALSE)
		 {
			alpha <- state$alpha
            m1 <- state$m1
            muclus <- state$muclus 
			ncluster <- state$ncluster
			psi1 <- state$psi1
			psiinv1 <- solve(psi1)
			k0 <- state$k0
			sigmaclus <- state$sigmaclus 
			ss <- state$ss
            if(nmissi==1) z <- state$z
		 }    

       #########################################################################################
       # working space
       #########################################################################################
         ccluster <- rep(0,nrec)
         cstrt <- matrix(0,nrow=nrec,ncol=nrec) 
         num <- matrix(0,nrow=npred,ncol=ngrid)
         denom <- rep(0,npred)  
         iflag <- rep(0,nvar) 
         muwork <- rep(0,nvar) 
         muwork2 <- rep(0,nvar) 
         prob <- rep(0,(nrec+100))
         s1 <- matrix(0,nvar,nvar)
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)
         sigmawork <- matrix(0,nrow=nvar,ncol=nvar)
         sigmawork2 <- matrix(0,nrow=nvar,ncol=nvar)
         sigworkinv <- matrix(0,nrow=nvar,ncol=nvar)
         theta <- rep(0,nvar)
         workm1 <- matrix(0,nrow=nvar,ncol=nvar)
         workm2 <- matrix(0,nrow=nvar,ncol=nvar)
         workm3 <- matrix(0,nrow=nvar,ncol=nvar)
         workmh1 <- rep(0,nvar*(nvar+1)/2) 
         workmh2 <- rep(0,nvar*(nvar+1)/2) 
         workv1 <- rep(0,nvar) 
         workv2 <- rep(0,nvar) 
         workv3 <- rep(0,nvar) 
		 ywork <- rep(0,nvar)

         iflagx <- rep(0,nx)
         workvx <- rep(0,nx) 
         workmx <- matrix(0,nrow=nx,ncol=nx) 

         fs <- rep(0,ngrid) 
         fm <- rep(0,npred) 

         worksam <- rep(0,nsave) 

         
         numcpo <- rep(0,nrec)
         denomcpo <- rep(0,nrec)
         
         nuvec <- c(nuvec,m1rand)

       #########################################################################################
       # calling the fortran code
       #########################################################################################


         foo <- .Fortran("dpdenregr",
				nrec       =as.integer(nrec),
                nx         =as.integer(nx),
				nvar       =as.integer(nvar),
				nmissi     =as.integer(nmissi),
				nmiss      =as.integer(nmiss),
				z          =as.double(z),
				missp      =as.integer(missp),
				npred      =as.integer(npred),
				xpred      =as.double(xpred),
				ngrid      =as.integer(ngrid),
				grid       =as.double(grid),
				a0b0       =as.double(a0b0),
				k0         =as.double(k0),
				nuvec      =as.integer(nuvec),
				s2inv      =as.double(s2inv),
				s2invm2    =as.double(s2invm2),
				psiinv2    =as.double(psiinv2),
				tau        =as.double(tau),
				mcmc       =as.integer(mcmcvec),
				nsave      =as.integer(nsave),
				cpo        =as.double(cpo),
				thetasave  =as.double(thetasave),
				denspm     =as.double(denspm),
				denspl     =as.double(denspl),
				densph     =as.double(densph),
				meanfpm    =as.double(meanfpm),
				meanfpl    =as.double(meanfpl),
				meanfph    =as.double(meanfph),
				alpha      =as.double(alpha),		
				m1         =as.double(m1),		
                muclus     =as.double(muclus),		 		
				ncluster   =as.integer(ncluster),
				psi1       =as.double(psi1),
				psiinv1    =as.double(psiinv1),
				s1         =as.double(s1),
				sigmaclus  =as.double(sigmaclus),
				ss         =as.integer(ss),
				ccluster   =as.integer(ccluster),
                cstrt      =as.integer(cstrt), 
				iflag      =as.integer(iflag),
				num        =as.double(num),
				denom      =as.double(denom),
				fs         =as.double(fs),
				fm         =as.double(fm),
				muwork     =as.double(muwork),
				prob       =as.double(prob),
				seed       =as.integer(seed),
				sigmawork  =as.double(sigmawork),
				sigworkinv =as.double(sigworkinv),
				theta      =as.double(theta),
				workm1     =as.double(workm1),
				workm2     =as.double(workm2),
				workm3     =as.double(workm3),
				workmh1    =as.double(workmh1),
				workmh2    =as.double(workmh2),
				workv1     =as.double(workv1),
				workv2     =as.double(workv2),
				workv3     =as.double(workv3),
				ywork      =as.double(ywork),
				iflagx     =as.integer(iflagx),
				workvx     =as.double(workvx),
				workmx     =as.double(workmx),
				worksam    =as.double(worksam),
				numcpo     =as.double(numcpo),
				denomcpo   =as.double(denomcpo),
				PACKAGE    ="DPpackage")

       #########################################################################################
       # save state
       #########################################################################################

         model.name <- "Bayesian Semiparametric Density Regression"		

         cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo<-cpom[,1]         
         fso<-cpom[,2]
       
         varnames <- colnames(as.matrix(x))
         if(is.null(varnames))
         {
            varnames<-all.vars(cl)[2]
         }
         varnames <- c("y",varnames)
	
         state <- list(
                  alpha=foo$alpha,
                  m1=matrix(foo$m1,nrow=nvar,ncol=1),
                  muclus=matrix(foo$muclus,nrow=nrec+100,ncol=nvar),
                  ncluster=foo$ncluster,
                  psi1=matrix(foo$psi1,nrow=nvar,ncol=nvar),
                  k0=foo$k0,
                  sigmaclus=matrix(foo$sigmaclus,nrow=nrec+100,ncol=nuniq),
                  ss=foo$ss,
                  z=matrix(foo$z,nrow=nrec,ncol=nvar) 
                  )

         pnames1 <- NULL
         for(i in 1:nvar)
		 {
			pnames1 <- c(pnames1,paste("m1",varnames[i],sep=":"))
		 }
         pnames2 <- "k0"
       
         pnames3 <- NULL
		 for(i in 1:nvar)
		 {
			for(j in i:nvar)
			{
				if(i==j)
				{
					tmp<-varnames[i]
				}
				else
				{
					tmp <- paste(varnames[i],varnames[j],sep="-")
				}   
				pnames3 <- c(pnames3,paste("psi1",tmp,sep=":"))
			}	
		 }

         pnames4 <- c("ncluster","alpha")
         pnames <- c(pnames1,pnames2,pnames3,pnames4)

         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)

         densp.m <- matrix(foo$denspm,nrow=npred,ncol=ngrid)
         densp.l <- matrix(foo$denspl,nrow=npred,ncol=ngrid)
         densp.h <- matrix(foo$densph,nrow=npred,ncol=ngrid)

         meanfp.m <- foo$meanfpm
         meanfp.l <- foo$meanfpl
         meanfp.h <- foo$meanfph

         colnames(thetasave) <- pnames 

         coeff <- apply(thetasave,2,mean)
 
         save.state <- list(thetasave=thetasave) 
 
         z <- list(call=cl,
                   y=y,
                   coefficients=coeff,
                   varnames=varnames,
                   modelname=model.name,
                   cpo=cpo,
                   fso=fso,  
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   save.state=save.state,
                   nrec=foo$nrec,
                   nvar=foo$nvar,
                   alpharand=alpharand,
                   psi1rand=psi1rand,
                   m1rand=m1rand,
                   k0rand=k0rand,
                   densp.m=densp.m,
                   densp.l=densp.l,
                   densp.h=densp.h,
                   meanfp.m=meanfp.m,
                   meanfp.l=meanfp.l,
                   meanfp.h=meanfp.h,
                   npred=npred,
                   ngrid=ngrid,
                   grid=grid,
                   compute.band=compute.band)
                 
         cat("\n\n")
         class(z)<-"DPcdensity"
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara, 2008
### Last modification: 25-05-2008.
###

"print.DPcdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Predictors:",x$nvar-1,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.DPcdensity"<-function(object, hpd=TRUE, ...) 
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

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo<-object$cpo

### Baseline Information

    mat<-NULL
    coef.p<-NULL
    
    dimen1<-object$nvar+object$nvar*(object$nvar+1)/2+1

    coef.p<-object$coefficients[1:dimen1]
    mat<-thetasave[,1:dimen1]
    
    if(dimen1>0){
    
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

    dimen1<-object$nvar+object$nvar*(object$nvar+1)/2+1
    
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
    ans$nvar<-object$nvar

    class(ans) <- "summaryDPcdensity"
    return(ans)
}


"print.summaryDPcdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Predictors:",x$nvar-1,"\n")        
    cat("\n\n")
    invisible(x)
}



"plot.DPcdensity"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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

   if(is(x, "DPcdensity"))
   {
      if(output=="density")
      {
      # Density estimation
		par(ask = ask)
		layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))

        for(i in 1:x$npred)
		{
               if(x$compute.band)
               {
                  title1 <- paste("Density Prediction #",i,sep=" ")           
                  plot(x$grid,x$densp.h[i,],main=title1,lty=2,type='l',lwd=2,xlab="y",ylab="density")
                  lines(x$grid,x$densp.l[i,],lty=2,lwd=2)
                  lines(x$grid,x$densp.m[i,],lty=1,lwd=3)
			   }
               else
               {
				  title1 <- paste("Density Prediction #",i,sep=" ")           
				  plot(x$grid,x$densp.m[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
			   }
		}

      }
      else
      {

        if(is.null(param))
        {
           pnames<-colnames(x$save.state$thetasave)
           n<-dim(x$save.state$thetasave)[2]
           cnames<-names(x$coefficients)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:(n-1))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                 fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
               fancydensplot1(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
        }
   
        else
        {

            pnames<-colnames(x$save.state$thetasave)
            n<-dim(x$save.state$thetasave)[2]
			poss<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }
            if (poss==0) 
			{
				stop("This parameter is not present in the original model.\n")
			}
	    
			par(ask = ask)
			layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
            title1<-paste("Trace of",pnames[poss],sep=" ")
            title2<-paste("Density of",pnames[poss],sep=" ")       
            plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
            if(pnames[poss]=="ncluster")
            {
                hist(x$save.state$thetasave[,poss],main=title2,xlab="values", ylab="probability",probability=TRUE)
            }
            else
            {
               fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            
        }
      }	
   }
}




