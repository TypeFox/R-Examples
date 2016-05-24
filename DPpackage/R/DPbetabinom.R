### DPbetabinom.R                   
### Fit a semiparametric beta binomial model using a DP prior.
###
### Copyright: Alejandro Jara, 2009 - 2012.
### Last modification: 09-03-2009.
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
###      Fernando Quintana
###      Departamento de Estadistica 
###      Facultad de Matematicas 
###      Pontificia Universidad Catolica de Chile 
###      Casilla 306, Correo 22
###      Santiago
###      Voice: +56-2-3544464  URL  : http://www.mat.puc.cl/~quintana
###      Fax  : +56-2-3547229  Email: quintana@mat.puc.cl
###


"DPbetabinom"<-
function(y,ngrid=100,prior,mcmc,state,status,data=sys.frame(sys.parent()),work.dir=NULL)
UseMethod("DPbetabinom")

DPbetabinom.default<-
function(y,
         ngrid=100,
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
         nrec <- nrow(y)
         q <- ncol(y)

         if(is.null(q))
         {
            stop("The respons must be a nrec*2 matrix.\n")      
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
         grid <- seq(from=0,to=1,length.out=ngrid)
         
       #########################################################################################
       # prior information
       #########################################################################################
         if(is.null(prior$a0))      
  	 {
            a0 <--1
  	    b0 <--1 
            alpha<-prior$alpha
         }
         else
         {
            a0 <- prior$a0
            b0 <- prior$b0
  	    alpha <- 1
  	 }

         a1 <- prior$a1
         b1 <- prior$b1

       #########################################################################################
       # mcmc specification
       #########################################################################################
         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave <- mcmc$nsave

       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2) 
         thetasave <- matrix(0,nrow=nsave,ncol=3)
         randsave <- matrix(0,nrow=nsave,ncol=(nrec+1))
         densm <- rep(0,ngrid)

       #########################################################################################
       # parameters depending on status
       #########################################################################################

         if(status==TRUE)
	 {
             ncluster <- 1
             ss <- rep(1,nrec)
             p <- rep(0,nrec+1)
             p[1] <- mean(y[,1]/y[,2])
   	 }
	 
      	 if(status==FALSE)
	 {
	    alpha <- state$alpha
            ncluster <- state$ncluster
	    ss <- state$ss
	    p <- state$p
	 }    

       #########################################################################################
       # working space
       #########################################################################################

         cstrt <- matrix(0,nrow=nrec,ncol=nrec)
         ccluster <- rep(0,nrec)
         prob <- rep(0,(nrec+1))
         lprob <- rep(0,(nrec+1))
         workcpo <- rep(0,nrec)

         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         y <- cbind(y[,1],y[,2]-y[,1])

         foo <- .Fortran("dpbetabinom",
  	 	nrec       =as.integer(nrec),
  	 	y          =as.double(y),
  	 	ngrid      =as.integer(ngrid),
  	 	grid       =as.double(grid),
  	 	a0         =as.double(a0),
  	 	b0         =as.double(b0),
  	 	a1         =as.double(a1),
  	 	b1         =as.double(b1),
 		ncluster   =as.integer(ncluster),
 		ss         =as.integer(ss),
		alpha      =as.double(alpha),		
                p          =as.double(p),
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
 		cpo        =as.double(cpo),
 		densm      =as.double(densm),
 		thetasave  =as.double(thetasave),
 		randsave   =as.double(randsave),
                ccluster   =as.integer(ccluster),
                cstrt      =as.integer(cstrt), 
                prob       =as.double(prob),
                lprob      =as.double(lprob),
                workcpo    =as.double(workcpo),
 		seed       =as.integer(seed),
		PACKAGE    ="DPpackage")

       #########################################################################################
       # save state
       #########################################################################################

         model.name <- "Bayesian Semiparametric Beta-Binomial Model"		

         cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo<-cpom[,1]         
         fso<-cpom[,2]
       
         state <- list(alpha=foo$alpha,
                       ncluster=foo$ncluster,
                       ss=foo$ss,
                       p=foo$p)

         pnames <- c("ncluster","alpha","LPML")
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=3)
         colnames(thetasave) <- pnames 

         coeff <- apply(thetasave,2,mean)[1:2]

         randsave<-matrix(foo$randsave,nrow=nsave,ncol=(nrec+1))

         densp.m <- foo$densm

         save.state <- list(thetasave=thetasave,
                            randsave=randsave) 
 
         z <- list(call=cl,
                   y=y,
                   coefficients=coeff,
                   modelname=model.name,
                   cpo=cpo,
                   fso=fso,  
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   save.state=save.state,
                   nrec=foo$nrec,
                   densp.m=densp.m,
                   ngrid=ngrid,
                   grid=grid)
                 
         cat("\n\n")
         class(z) <- "DPbetabinom"
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara, 2009
### Last modification: 09-03-2009.
###

"print.DPbetabinom"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\n\n")
    invisible(x)
}


"summary.DPbetabinom"<-function(object, hpd=TRUE, ...) 
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

### Precision parameter

    if(is.null(object$prior$a0))
    {
      coef.p <- object$coefficients[1]
      mat <- matrix(thetasave[,1],ncol=1)
    }
    else
    {
      coef.p <- object$coefficients
      mat <- thetasave[,1:2]

    }  

    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

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

    ans$prec <- coef.table
    ans$nrec <- object$nrec

    class(ans) <- "summaryDPbetabinom"
    return(ans)
}


"print.summaryDPbetabinom"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$prec)) {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}



"plot.DPbetabinom"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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

   if(is(x, "DPbetabinom"))
   {
      if(output=="density")
      {
      # Density estimation
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))

        title1 <- "Posterior density estimate"
        plot(x$grid,x$densp.m,lwd=2,type="l",lty=1,main=title1,xlab="values",ylab="density")

      }
      else
      {

        if(is.null(param))
        {
           pnames <- colnames(x$save.state$thetasave)
           cnames <- names(x$coefficients)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:1)
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
               title1<-paste("Trace of",pnames[2],sep=" ")
               title2<-paste("Density of",pnames[2],sep=" ")       
               plot(x$save.state$thetasave[,2],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,2],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-paste("Trace Plot of",pnames[3],sep=" ")
           plot(x$save.state$thetasave[,3],type='l',main=title1,xlab="MCMC scan",ylab="LPML")
        }
   
        else
        {

            pnames <- colnames(x$save.state$thetasave)
            n <- ncol(x$save.state$thetasave)
	    poss <- 0 
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
            title1 <- paste("Trace of",pnames[poss],sep=" ")
            title2 <- paste("Density of",pnames[poss],sep=" ")       
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




