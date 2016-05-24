### LDBDPdensity.R                    
### Fit a linear dependent Bernstein-Dirichlet model for conditional density estimation.
###
### Copyright: Felipe Barrientos and Alejandro Jara, 2010-2012.
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
### The authors' contact information:
###
###      Felipe Barrientos
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Email: afbarrie@mat.puc.cl
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
###

"LDBDPdensity"<-
function(formula,xpred,prior,mcmc,state,status,ngrid=100,grid=NULL,compute.band=FALSE,type.band="PD",data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("LDBDPdensity")

"LDBDPdensity.default"<-
function(formula,
         xpred,
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
		 y <- model.response(mf,"numeric")
		 nrec <- length(y)
		 x <- model.matrix(formula)
		 p <- ncol(x)

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

         npred <- nrow(xpred)

         if(is.null(grid))
         {
			grid <- seq(0,1,len=ngrid)
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
       # Prior information
       #########################################################################################

		 if(is.null(prior$maxn))
		 {
			maxn <- 50
		 }
		 else
		 {
		    maxn <- prior$maxn
		 }

		 alpha  <- prior$alpha


		 s0inv <- solve(prior$S0)
		 s0invm <- s0inv%*%prior$m0

		 if(is.null(prior$lambda))
	     {
			if(is.null(prior$k))
			{
				cat("\n *** Error: k MUST be specified.\n")
				return(-1)
			}
            kk <- prior$k
            lambda <- -1
		 }
		 else
		 {
			lambda <- prior$lambda
			kk <- lambda
		 }

		 if(is.null(prior$nu))
		 {
			if(is.null(prior$tau1))
			{
				cat("\n *** Error: tau1 and tau2 MUST be specified.\n")
				return(-1)
			}
			if(prior$tau1<0)
			{
				cat("\n *** Error: tau1 and tau2 MUST be positive.\n")
				return(-1)
			}
			tau1   <- prior$tau1
			tau2   <- prior$tau2
			gp <- 2*nrec
			
			nu <- -3
			psiinv <- diag(1,p)
		 }
		 else
		 {
			tau1 <- -1
			tau2 <- -1
			gp <- 1

			nu <- prior$nu
			psiinv <- prior$psiinv 
		 }

       #########################################################################################
       # mcmc specification
       #########################################################################################

         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,cband,tband)
         nsave <- mcmc$nsave

         if(is.null(mcmc$slicebeta))
         {
			slicebeta <- 0.5
		 }
		 else
		 {
			slicebeta <- mcmc$slicebeta
		 }

         if(is.null(mcmc$sliceg))
         {
			slicegp <- 5
		 }
		 else
		 {
			slicegp <- mcmc$sliceg
		 }

         if(is.null(mcmc$slicev))
         {
			slicev <- 0.01
		 }
		 else
		 {
			slicev <- mcmc$slicev
		 }

		 slice  <- c(slicebeta,slicegp,slicev) 

       #########################################################################################
       # output
       #########################################################################################
		 acrate <- c(0,0)
		 cpo <- matrix(0,nrow=nrec,ncol=2)
         denspm <- matrix(0,nrow=npred,ncol=ngrid)
         denspl <- matrix(0,nrow=npred,ncol=ngrid)
         densph <- matrix(0,nrow=npred,ncol=ngrid)

		 meanfpm <- rep(0,npred)
		 meanfpl <- rep(0,npred)
		 meanfph <- rep(0,npred)


         thetasave <- matrix(0,nrow=nsave,ncol=2+p+p*(p+1)/2)
		 randsave <- matrix(0,nrow=nsave,ncol=p*maxn+maxn)

       #########################################################################################
       # parameters depending on status
       #########################################################################################

         if(status==TRUE)
         {
           beta <- matrix(0,nrow=maxn,ncol=p)
           mub <- solve(t(x)%*%x)%*%t(x)%*%y 
           sb <- solve(t(x)%*%x)
           v <- c(rep(0.5,maxn-1),1)
         }

         if(status==FALSE) 
         {
           if(!is.null(prior$lambda)){ kk <- state$k }
		   if(is.null(prior$nu)){ gp <- state$g }
           beta <- state$beta
           mub <- state$mub
           sb <- state$Sb
           v <- state$v
         }    

       #########################################################################################
       # working space
       #########################################################################################
         
		 seed <- c(sample(1:29000,1),sample(1:29000,1))

		 iflag <- rep(0,p)
		 sbinv <- matrix(0,nrow=p,ncol=p)
		 workm1 <- matrix(0,nrow=p,ncol=p)
		 workv1 <- rep(0,p)
		 workv2 <- rep(0,p)
		 workmh1 <- rep(0,p*(p+1)/2)
		 workmh2 <- rep(0,p*(p+1)/2)

		 betal <- matrix(0,nrow=maxn,ncol=p)
		 betar <- matrix(0,nrow=maxn,ncol=p)
	 	 beta1 <- matrix(0,nrow=maxn,ncol=p)
		 vl <- rep(0,maxn)
		 vr <- rep(0,maxn)
		 v1 <- rep(0,maxn) 

		 workdpw <- rep(0,maxn+1)
		 weight  <- rep(0,maxn)

		 fw <- rep(0,nrec)
		 fw2 <- matrix(0,nrow=npred,ncol=ngrid)

		 fm <- rep(0,npred)

		 fs <- rep(0,ngrid)

		 worksam <- rep(0,nsave)


       #########################################################################################
       # calling the fortran code
       #########################################################################################

          foo <- .Fortran("ldbdpdensity",
                 y			=as.double(y),
                 x			=as.double(x),
                 nrec		=as.integer(nrec),
                 p			=as.integer(p),
                 npred		=as.integer(npred),
                 ngrid		=as.integer(ngrid),
                 grid		=as.double(grid),
                 xpred		=as.double(xpred),
                 maxn		=as.integer(maxn),
                 nu			=as.integer(nu),
                 alpha		=as.double(alpha),
                 lambda		=as.double(lambda),
                 tau1		=as.double(tau1),
                 tau2		=as.double(tau2),
                 psiinv		=as.double(psiinv),
                 s0invm		=as.double(s0invm),
                 s0inv		=as.double(s0inv),
                 kk			=as.integer(kk),
                 gp			=as.double(gp),
                 beta		=as.double(beta),
                 mub		=as.double(mub),
                 sb			=as.double(sb),
                 v			=as.double(v),
                 mcmc		=as.integer(mcmcvec),
                 nsave		=as.integer(nsave),
                 slice		=as.double(slice),
                 acrate		=as.double(acrate),
                 thetasave	=as.double(thetasave),
				 randsave	=as.double(randsave),
                 fmean		=as.double(denspm),
                 flow		=as.double(denspl),
                 fupp		=as.double(densph),
				 meanfpm   = as.double(meanfpm),
				 meanfpl   = as.double(meanfpl),
				 meanfph   = as.double(meanfph),
				 cpo        =as.double(cpo),
                 seed		=as.integer(seed),
                 iflag		=as.integer(iflag),
                 sbinv		=as.double(sbinv),
                 workm1		=as.double(workm1),
                 workv1		=as.double(workv1),
                 workv2		=as.double(workv2),
                 workmh1	=as.double(workmh1),
                 workmh2	=as.double(workmh2),
                 betal		=as.double(betal),
                 betar		=as.double(betar),
                 beta1		=as.double(beta1),
                 vl			=as.double(vl),
                 vr			=as.double(vr),
                 v1			=as.double(v1),
                 workdpw	=as.double(workdpw),
                 weight		=as.double(weight),
                 fw			=as.double(fw),
                 fw2		=as.double(fw2),
                 fs			=as.double(fs),
				 fm        = as.double(fm),
				 worksam	=as.double(worksam),
                 PACKAGE="DPpackage")
         
        
       #########################################################################################
       # save state
       #########################################################################################

         if(!is.null(work.dir))
         {
            cat("\n\n Changing working directory back to ",old.dir,"\n")
            setwd(old.dir)
         }

         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo <- cpom[,1]         
         fso <- cpom[,2]

         model.name <-"Bounded density regression using dependent Bernstein polynomials"

		 state <- list(k=foo$kk,
					   g=foo$gp,
					   beta=matrix(foo$beta,nrow=maxn,ncol=p),
					   mub=foo$mub,
					   Sb=matrix(foo$sb,nrow=p,ncol=p),
					   v=foo$v)


         denspm <- matrix(foo$fmean,nrow=npred,ncol=ngrid)
         denspl <- matrix(foo$flow,nrow=npred,ncol=ngrid)
         densph <- matrix(foo$fupp,nrow=npred,ncol=ngrid)

		 meanfpm <- foo$meanfpm
		 meanfpl <- foo$meanfpl
		 meanfph <- foo$meanfph

		 thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=2+p+p*(p+1)/2)
		 randsave <- matrix(foo$randsave,nrow=nsave,ncol=maxn*p+maxn)

		 coeffname <- dimnames(x)[[2]]

		 pnames<-vector()
		 rnames<-vector()


		 count <- 0
		 count <- count+1
		 pnames[count] <- "k"

		 for(i in 1:p)
		 {
			count <- count+1
			pnames[count] <- paste("mub",coeffname[i],sep = "-")
		 }

		 for(i in 1:p)
		 {
            for(j in i:p)
            {
				count<-count+1
				pnames[count] <- paste("Sb",coeffname[i],coeffname[j],sep = "-")
			}
		 }

		 count <- count+1
		 pnames[count] <- "g"


		 count <- 0	
		 for(i in 1:maxn)
		 {
			 for(j in 1:p)
             {
				 count<-count+1
				 rnames[count]<-paste(coeffname[j],i,sep="-")
			 }
		 }
		 for(j in 1:maxn)
		 {
			count<-count+1
			rnames[count]<-paste("v",j,sep = "")
		 }


		 colnames(thetasave) <- pnames 
		 colnames(randsave) <- rnames 

		 coeff <- apply(thetasave,2,mean)

         nn <- length(coeff)
		 if(!is.null(prior$nu)) coeff <- coeff[-nn]
         
         save.state <- list(thetasave=thetasave,
							randsave=randsave)

		 z <- list(	slice=foo$acrate[1],
					acrate=foo$acrate[2],
					modelname=model.name,
					call=cl,
					cpo=cpo,
					coefficients=coeff,
					fso=fso,
					prior=prior,
					mcmc=mcmc,
					nrec=foo$nrec,
					p=foo$p,
					x=x,
					tau1=tau1,
					ngrid=ngrid,
					npred=npred,
					xpred=xpred,
					grid=grid,
					densp.m=denspm,
					densp.l=denspl,
					densp.h=densph,
					meanfp.m=meanfpm,
					meanfp.l=meanfpl,
					meanfp.h=meanfph,
					state=state,
					save.state=save.state,
					work.dir=work.dir,
                    compute.band=compute.band)

		 cat("\n\n")
		 class(z)<-c("LDBDPdensity")
		 z 
}



###                    
### Tools
###
### Copyright: Alejandro Jara, 2008
### Last modification: 02-06-2008.
###

"print.LDBDPdensity" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors:",x$p,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.LDBDPdensity"<-function(object, hpd=TRUE, ...) 
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
         alpha <- 0.05
         vec <- x
         n <- length(x)         
         alow <- rep(0,2)
         aupp <- rep(0,2)
         a <- .Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

    thetasave <- object$save.state$thetasave

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo <- object$cpo

### Bernstein polynomial degree

    dimen1 <- 1
    coef.p <- object$coefficients[dimen1]
    mat <- matrix(thetasave[,dimen1],ncol=1)

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

	ans$bernstein <- coef.table

### Baseline Information

    mat <- NULL
    coef.p <- NULL
	p <- object$p

    dimen2 <- p+p*(p+1)/2
    
    coef.p <- object$coefficients[(dimen1+1):(dimen1+dimen2)]
    mat <- thetasave[,(dimen1+1):(dimen1+dimen2)]

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

	ans$base <- coef.table

### g-prior
 
    if(!is.null(object$prior$nu))
    {
		ans$gprior <- NULL
    }
    else
    {
		coef.p <- object$coefficients[(dimen1+dimen2+1)]
		mat <- matrix(thetasave[,(dimen1+dimen2+1)],ncol=1)
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
		ans$gprior <- coef.table

	}
    ans$acrate <- object$acrate
    ans$nrec <- object$nrec
    ans$p <- object$p

    class(ans) <- "summaryLDBDPdensity"
    return(ans)
}


"print.summaryLDBDPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 

    if (length(x$bernstein)) {
        cat("\nDegree of Bernstein polynomial:\n")
        print.default(format(x$bernstein, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No degree of Bernstein polynomial\n")


    if (length(x$base)) {
        cat("\nDP baseline distribution:\n")
        print.default(format(x$base, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No baseline parameters\n")

    if (length(x$gprior)) {
        cat("\nG-prior:\n")
        print.default(format(x$gprior, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

	cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors:",x$p,"\n")        
    cat("\n\n")
    invisible(x)
}





"plot.LDBDPdensity"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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
		alpha <- 0.05
		alow <- rep(0,2)
		aupp <- rep(0,2)
		n <- length(x)
		a <- .Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(x),
		                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
		xlinf <- a$alow[1]            
		xlsup <- a$aupp[1]            
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

        
	xpol <- c(xlinf,xlinf,densx[densx>=xlinf & densx <=xlsup],xlsup,xlsup)
	ypol <- c(0,ylinf,densy[densx>=xlinf & densx <=xlsup] ,ylsup,0)
             
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


   if(is(x, "LDBDPdensity"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:(n-1))
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="k")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }
		   if(is.null(x$prior$nu))
		   {
			   title1 <- paste("Trace of",pnames[n],sep=" ")
               title2 <- paste("Density of",pnames[n],sep=" ")       
               plot(ts(x$save.state$thetasave[,n]),main=title1,xlab="MCMC scan",ylab=" ")
			   fancydensplot1(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
		   }

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
        }
   }

}





