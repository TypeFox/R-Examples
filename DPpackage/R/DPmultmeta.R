### DPmultmeta.R                   
### Fit a semiparametric random effects model for multivariate meta-analysis 
### using a Dirichlet Process prior for the distribution of 
### the random effects.
###
### Copyright: Alejandro Jara and Peter Mueller, 2008-2012.
###
### Last modification: 02-06-2008.
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
###      Peter Mueller
###      Department of Mathematics
###      The University of Texas Austin
###      1, University Station, C1200 
###      Austin TX 78712, USA
###      Voice: (512) 471-7168  URL  : http://www.math.utexas.edu/users/pmueller
###      Fax  : (512) 471-9038  Email: pmueller@math.utexas.edu
###

"DPmultmeta"<-
function(y,asymvar,prior,mcmc,state,status,data=sys.frame(sys.parent()))
UseMethod("DPmultmeta")

"DPmultmeta.default"<-
function(y,
         asymvar, 
         prior,
         mcmc,
         state,
         status, 
         data=sys.frame(sys.parent()))
{
       #########################################################################################
       # call parameters
       #########################################################################################
	 cl <- match.call()

       #########################################################################################
       # data structure
       #########################################################################################

         nameresp <- colnames(y)
         nvar <- dim(y)[2]
         nrec <- dim(y)[1]

         if(nvar < 2)
         {
            stop("Use the function DPMmeta for univariate meta-analysis")
         }   

  	 sigma2e <- asymvar
         nrecs <- dim(sigma2e)[1]
         nvars <- dim(sigma2e)[2]

         if(nrec != nrecs)
         {
            stop("Different number of subjects in the response and the assymptotic variance matrix")
         }   

         nuniq <- nvar*(nvar+1)/2
         if(nuniq != nvars)
         {
            stop("Different dimension in the response vector and the corresponding assymptotic variance")
         }   

       #########################################################################################
       # prior information
       #########################################################################################

  	 if(is.null(prior$a0))
  	 {
  	    a0 <--1
  	    b0 <--1 
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

  	 if(is.null(prior$m2))
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
       		  m1[i] <- mean(y[,i])+rnorm(1,0,100)
       	      }
              m1rand <- 1
         }     

  	 if(is.null(prior$nu))
  	 {
              s1 <- matrix(prior$s1,nvar,nvar)
              psiinv <- matrix(0,nvar,nvar)
              s1rand <- 0
              nu <- -1
  	 }
  	 else
  	 {
  	      s1 <- matrix(var(y),nvar,nvar)
              psiinv <- matrix(prior$psiinv,nvar,nvar)
  	      s1rand <- 1
              nu <- prior$nu
  	 }
  	 

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
            nsave <-mcmc$nsave
         }

       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2)
         thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+2)
         randsave <- matrix(0,nrow=nsave,ncol=(nrec+1)*nvar)

       #########################################################################################
       # parameters depending on status
       #########################################################################################
         nuniq <- nvar*(nvar+1)/2
         
    	 if(status==TRUE)
	 {
            muclus <- matrix(0,nrow=nrec+1,ncol=nvar)
            for(i in 1:1)
            {
                for(j in 1:nvar)
                {
		    muclus[i,j] <- m1[j]
                }
            }
                ncluster <- 1
                ss <- rep(1,nrec)
   	 }
	 
      	 if(status==FALSE)
	 {
	        alpha <- state$alpha
                m1 <- state$m1
                s1 <- state$s1
                muclus <- state$muclus 
	        ncluster <- state$ncluster
	        ss <- state$ss
	 }    
         
       #########################################################################################
       # working space
       #########################################################################################
         ccluster <- rep(0,nrec)
         cstrt <- matrix(0,nrow=nrec,ncol=nrec) 
         iflag <- rep(0,nvar) 
         prob <- rep(0,(nrec+1))
         sigma2ei <- matrix(0,nrow=nrec,ncol=nvar*(nvar+1)/2) 
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)
         theta <- rep(0,nvar)
         s1inv <- matrix(0,nrow=nvar,ncol=nvar)
         s1invm1 <- rep(0,nvar)
         workm1 <- matrix(0,nrow=nvar,ncol=nvar)
         workm2 <- matrix(0,nrow=nvar,ncol=nvar)
         workm3 <- matrix(0,nrow=nvar,ncol=nvar)
         workmh1 <- rep(0,nvar*(nvar+1)/2)
         workmh2 <- rep(0,nvar*(nvar+1)/2)
         workv1 <- rep(0,nvar) 
         workv2 <- rep(0,nvar) 
	 ywork <- rep(0,nvar)

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         foo <- .Fortran("dpmultmeta",
                nrec       =as.integer(nrec),
                nvar       =as.integer(nvar),
                y          =as.double(y),
                sigma2e    =as.double(sigma2e),	 
                a0b0       =as.double(a0b0),
                m1rand     =as.integer(m1rand),
                s2inv      =as.double(s2inv),
                s2invm2    =as.double(s2invm2),
                nu         =as.integer(nu),
                s1rand     =as.integer(s1rand),
                psiinv     =as.double(psiinv),	  		
                mcmc       =as.integer(mcmcvec),
                nsave      =as.integer(nsave),
                cpo        =as.double(cpo),
                randsave   =as.double(randsave),
                thetasave  =as.double(thetasave),
                alpha      =as.double(alpha),
                m1         =as.double(m1),
                s1         =as.double(s1),
                ncluster   =as.integer(ncluster),
                muclus     =as.double(muclus),
                ss         =as.integer(ss),
                ccluster   =as.integer(ccluster),
                cstrt      =as.integer(cstrt),
                iflag      =as.integer(iflag),
                prob       =as.double(prob),
                sigma2ei   =as.double(sigma2ei),
                seed       =as.integer(seed),
                s1inv      =as.double(s1inv),
                s1invm1    =as.double(s1invm1),
                theta      =as.double(theta),
                workm1     =as.double(workm1),
                workm2     =as.double(workm2),
                workm3     =as.double(workm3),
                workmh1    =as.double(workmh1),
                workmh2    =as.double(workmh2),
                workv1     =as.double(workv1),
                workv2     =as.double(workv2),
                ywork      =as.double(ywork),
		PACKAGE    ="DPpackage")	


       #########################################################################################
       # save state
       #########################################################################################

	 model.name<-"Bayesian semiparametric random effects model for multivariate meta-analysis"

         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+2)
         randsave <- matrix(foo$randsave,nrow=nsave,ncol=(nrec+1)*nvar)
         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo <- cpom[,1]         
         fso <- cpom[,2]

	 varnames <- colnames(y)
         if(is.null(varnames))
         {
               varnames <- all.vars(cl)[1:nvar]
         }
         indip <- rep(0,(nvar+nvar*(nvar+1)/2+2))         

         coeff <- NULL
         pnames1 <- NULL
         pnames2 <- NULL
         renames <- NULL  

         for(i in 1:nvar)
         {
             pnames2 <- c(pnames2,paste("m1",varnames[i],sep="-"))
         }
         for(i in 1:nvar)
         {
             for(j in i:nvar)
             {
                 tmp <- paste(varnames[i],varnames[j],sep=":")
                 pnames2 <- c(pnames2,paste("s1",tmp,sep="-"))
             }
         }
         pnames2 <- c(pnames2,"ncluster","alpha")

         for(i in 1:nrec)
         {
             for(j in 1:nvar)
             {
                 tmp1 <- paste("m1",varnames[i],sep="-") 
                 tmp2 <- paste("Subject=",i,sep="")
                 renames <- c(renames,paste(tmp2,tmp1,sep=";"))
             } 
         } 

         for(j in 1:nvar)
         {
             tmp1 <- paste("m1",varnames[i],sep="-") 
             renames <- c(renames,paste("Prediction",tmp1,sep=";"))
         } 

         colnames(thetasave) <- pnames2
         colnames(randsave) <- renames                   

         count <- 0
         if(m1rand==1)
         {
            for(i in 1:nvar)
	    {
                count <- count + 1
	        coeff <- c(coeff,mean(thetasave[,i]))
	        pnames1 <- c(pnames1,paste("m1",varnames[i],sep="-"))
	        indip[i] <- 1
            }
         }

         count <- nvar
         if(s1rand==1)
         {
            for(i in 1:nvar)
	    {
                for(j in i:nvar)
                {
                    count <- count + 1
	            coeff <- c(coeff,mean(thetasave[,count]))
                    tmp <- paste(varnames[i],varnames[j],sep=":")
   	            pnames1 <- c(pnames1,paste("s1",tmp,sep="-"))
	            indip[count] <- 1
                }
            }
         }

         count <- nvar+nvar*(nvar+1)/2 + 1       
         coeff <- c(coeff,mean(thetasave[,count]))         
         pnames1 <- c(pnames1,"ncluster")
         indip[count] <- 1

         count <- nvar+nvar*(nvar+1)/2 + 1
         if(alpharand==1)
         {
            count <- count + 1
            coeff <- c(coeff,mean(thetasave[,count])) 
            pnames1 <- c(pnames1,"alpha")
            indip[count] <- 1
         }

         names(coeff) <- pnames1
 
         save.state <- list(thetasave=thetasave,
                            randsave=randsave)


         state <- list(alpha=foo$alpha,
                       m1=matrix(foo$m1,nrow=nvar,ncol=1),
                       s1=matrix(foo$s1,nrow=nvar,ncol=nvar),
                       muclus=matrix(foo$muclus,nrow=nrec+1,ncol=nvar),
                       ncluster=foo$ncluster,
                       ss=foo$ss)
         

	  z<-list(call=cl,
                  y=y,
                  asymvar=asymvar,
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
                  s1rand=s1rand,
                  m1rand=m1rand,
                  coefficients=coeff,
                  indip=indip) 
         cat("\n\n")        

         class(z)<-c("DPmultmeta")
         return(z) 
}


###                    
### Tools: print, summary, plot
###
### Copyright: Alejandro Jara, 2008
### Last modification: 02-06-2008.


"print.DPmultmeta"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Variables:",x$nvar,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.DPmultmeta"<-function(object, hpd=TRUE, ...) 
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
    
    dimen1 <- object$nvar+object$nvar*(object$nvar+1)/2
    
    for(i in 1:dimen1)
    {
       if(object$indip[i]==1)
       {
           coef.p <- c(coef.p,object$coefficients[i])
           mat <- cbind(mat,thetasave[,i])
       }   
    }

    dimen1 <- dim(mat)[2] 

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

    dimen1<-object$nvar+object$nvar*(object$nvar+1)/2
    
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

    class(ans) <- "summaryDPmultmeta"
    return(ans)
}


"print.summaryDPmultmeta"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Variables:",x$nvar,"\n")        
    cat("\n\n")
    invisible(x)
}


"plot.DPmultmeta"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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

"bivk"<-function(x, y, h, n = 25, lims = c(range(x), range(y))) 
{
    nx <- length(x)
    if (length(y) != nx) 
        stop("Data vectors must be the same length")
    gx <- seq(lims[1], lims[2], length = n)
    gy <- seq(lims[3], lims[4], length = n)
    if (missing(h)) 
        h <- c(band(x), band(y))
    h <- h/4
    ax <- outer(gx, x, "-")/h[1]
    ay <- outer(gy, y, "-")/h[2]
    z <- matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n, 
        nx))/(nx * h[1] * h[2])
    return(list(x = gx, y = gy, z = z))
}

"band"<-function(x) 
{
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2] - r[1])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}

   if(is(x, "DPmultmeta"))
   {

      if(output=="density")
      {

      # Density estimation
	
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))
	start <- x$nrec*x$nvar

        for(i in 1:x$nvar)
 	{
            tmp <- paste("m1",x$varnames[i],sep="-")
	    title1 <- paste("Density of",tmp,sep=' ')
            plot(density(x$save.state$randsave[,(start+i)]),type="l",lwd=2,,xlab="values", ylab="density",main=title1)
	}
	
        for(i in 1:(x$nvar-1))
	{
            for(j in (i+1):x$nvar)
            {
 	        varsn1 <- paste("m1",x$varnames[i],sep="-")
 	        varsn2 <- paste("m1",x$varnames[j],sep="-")
                tmp <- paste(varsn1,varsn2,sep=":")
	        title1<-paste("Density of ",tmp,sep='')
                xx<-x$save.state$randsave[,(start+i)]	    
                yy<-x$save.state$randsave[,(start+j)]	    
	        est<-bivk(xx,yy,n=200)
	        contour(est,main=title1,xlab=varsn1,ylab=varsn2)
	        persp(est,theta=-30,phi=15,expand = 0.9, ltheta = 120,main=title1,
	                 xlab=varsn1,ylab=varsn2,zlab="density")
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
               if(x$indip[i]==1)
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
           }
           
           if(!is.null(x$prior$a0))
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


