### TDPdensity.R                   
### For density estimation using a Triangular-Dirichlet prior
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 05-07-2007.
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

TDPdensity<-function(y,support=3,transform=1,ngrid=1000,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("TDPdensity")

TDPdensity.default<-function(y,support=3,transform=1,ngrid=1000,prior,mcmc,state,status,data,na.action=na.fail)
{
         #########################################################################################
         # call parameters
         #########################################################################################
           cl <- match.call()
           resp<-na.action(as.matrix(y))	
           varnames<-all.vars(cl)[1]
	  
         #########################################################################################
         # data structure
         #########################################################################################
           nrec<-dim(resp)[1]
           nvar<-dim(resp)[2]
           
           if(nvar>1)
           {
             stop("This function can only be used for univariate density estimation.\n")      
           }
           
           resp<-as.vector(resp)
           grid<-seq(0,1,length=ngrid)

           if(transform>1 && support==2)
           {
             stop("The parametric transformation can only be used for the real line support.\n")      
           }

           if(transform>1 && support==1)
           {
             stop("The parametric transformation can only be used for the real line support.\n")      
           }
           
           if(transform==1)
           {
              if(support==3)
              {
                 left<-min(resp)-0.5*sqrt(var(resp))
                 right<-max(resp)+0.5*sqrt(var(resp))
                 x<-(resp-left)/(right-left)
                 jacob<-1.0/(right-left)
                 grids<- left+grid*(right-left)
              }
              if(support==1)
              {
                 left<-0
                 right<-1
                 x<-resp
                 jacob<-1
                 grids<-grid
              }
              if(support==2)
              {
                 left<-0
                 right<-max(resp)+0.5*sqrt(var(resp))
                 x<-(resp-left)/(right-left)
                 jacob<-1.0/right
                 grids<- grid*right
              }
           }
           else
           {
              left<-min(resp)-0.5*sqrt(var(resp))
              right<-max(resp)+0.5*sqrt(var(resp))
              grids<- left+grid*(right-left)
              x<-rep(0,nrec)
              mu<-mean(resp)
              sigma2<-var(resp)
           }
           
           
         #########################################################################################
         # prior information
         #########################################################################################

  	   if(is.null(prior$aa0))
  	   {
  	      aa0<--1
  	      ab0<--1 
  	      alpha<-prior$alpha
  	      alpharand<-0
  	   }
           else
           {
              aa0<-prior$aa0
  	      ab0<-prior$ab0
  	      alpha<-1
  	      alpharand<-1
  	   }
  	   a0b0<-c(aa0,ab0)

           kmax<-prior$kmax
           a0<-prior$a0
           b0<-prior$b0
         

         #########################################################################################
         # mcmc specification
         #########################################################################################
          if(missing(mcmc))
          {
             nburn <- 1000
             nsave <- 1000
             nskip <- 0
             ndisplay <- 100
             mcmcvec<-c(nburn,nskip,ndisplay)
         }
         else
         {
             mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
             nsave<-mcmc$nsave
         }    

         #########################################################################################
         # output
         #########################################################################################

           fun<-rep(0,ngrid)
           cpo<-matrix(0,nrow=nrec,ncol=2)
           thetasave<-matrix(0,nrow=nsave,ncol=3)
           randsave<-matrix(0,nrow=nsave,ncol=(nrec+1))
         
         #########################################################################################
         # parameters depending on status
         #########################################################################################
         
    	   if(status==TRUE)
	   {
                  yclus<-rbeta(nrec,a0,b0)  
                  k<-kmax
                  ncluster<-nrec
                  ss<-seq(1,nrec)
                  if(!is.null(prior$mu))mu<-prior$mu
                  if(!is.null(prior$sigma2))sigma2<-prior$sigma2
   	   }
	 
      	   if(status==FALSE)
	   {
	          alpha<-state$alpha
                  k<-state$k
	          ncluster<-state$ncluster
	          yclus<-state$yclus
	          ss<-state$ss
	   }    

         #########################################################################################
         # working space
         #########################################################################################
           cstrt<-matrix(0,nrow=nrec,ncol=nrec)
           ccluster<-rep(0,nrec)
           prob<-rep(0,nrec+1)
           probk<-rep(0,kmax+1)
           seed<-c(sample(1:29000,1),sample(1:29000,1))
           y<-rep(0,nrec)

         #########################################################################################
         # calling the fortran code
         #########################################################################################

           if(transform==1)
           {
              foo <- .Fortran("tdpdensity",
                    nrec       =as.integer(nrec),
                    jacob      =as.double(jacob),  	 	
                    x          =as.double(x),
                    a0b0       =as.double(a0b0),
                    a0         =as.double(a0),
                    b0         =as.double(b0),  	 	
                    kmax       =as.integer(kmax),
                    k          =as.integer(k),
                    ncluster   =as.integer(ncluster),
                    ss         =as.integer(ss),
                    alpha      =as.double(alpha),
                    yclus      =as.double(yclus),
                    mcmc       =as.integer(mcmcvec),
                    nsave      =as.integer(nsave),
                    cpo        =as.double(cpo),
                    randsave   =as.double(randsave),
                    thetasave  =as.double(thetasave), 		
                    ngrid      =as.integer(ngrid),
                    grid       =as.double(grid),
                    fun        =as.double(fun),
                    seed       =as.integer(seed),
                    cstrt      =as.integer(cstrt), 		
                    ccluster   =as.integer(ccluster),
                    prob       =as.double(prob),
                    probk      =as.double(probk),
                    y          =as.double(y),
                    PACKAGE    ="DPpackage")
           }
           else
           {
              
              typet<-transform-1

              foo <- .Fortran("tdpdensitypl",
                    nrec       =as.integer(nrec),
                    resp       =as.double(resp),
                    a0b0       =as.double(a0b0),
                    a0         =as.double(a0),
                    b0         =as.double(b0),  	 	
                    kmax       =as.integer(kmax),
                    typet      =as.integer(typet),
                    mu         =as.double(mu),
                    sigma2     =as.double(sigma2),
                    k          =as.integer(k),
                    ncluster   =as.integer(ncluster),
                    ss         =as.integer(ss),
                    alpha      =as.double(alpha),
                    yclus      =as.double(yclus),
                    mcmc       =as.integer(mcmcvec),
                    nsave      =as.integer(nsave),
                    cpo        =as.double(cpo),
                    randsave   =as.double(randsave),
                    thetasave  =as.double(thetasave), 		
                    ngrid      =as.integer(ngrid),
                    grid       =as.double(grids),
                    fun        =as.double(fun),
                    seed       =as.integer(seed),
                    cstrt      =as.integer(cstrt), 		
                    ccluster   =as.integer(ccluster),
                    prob       =as.double(prob),
                    probk      =as.double(probk),
                    x          =as.double(x),
                    y          =as.double(y),
                    PACKAGE    ="DPpackage")
           }


         #########################################################################################
         # save state
         #########################################################################################

           cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)
           cpo<-cpom[,1]         
           fso<-cpom[,2]

           model.name<-"Bayesian semiparametric density estimation"		
         
           state <- list(alpha=foo$alpha,
                         yclus=foo$yclus, 
                         ncluster=foo$ncluster,
                         ss=foo$ss)
                       
           randsave<-matrix(foo$randsave,nrow=nsave,ncol=(nrec+1))
           thetasave<-matrix(foo$thetasave,nrow=nsave,ncol=3)
 
           pnames<-c("k","ncluster","alpha")
           colnames(thetasave)<-pnames
           coeff<-apply(thetasave,2,mean)

           save.state <- list(thetasave=thetasave,randsave=randsave)
           
  	   z<-list(call=cl,
  	           y=resp,
  	           x=x,
  	           varnames=varnames,
  	           cpo=cpo,
  	           fso=fso,
  	           modelname=model.name,
  	           coefficients=coeff,
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   save.state=save.state,
                   nrec=foo$nrec,
                   grid=grids,
                   grids=grid,
                   fun=foo$fun)
                 
          cat("\n\n")
 	  class(z)<-"TDPdensity"
  	  return(z)
}



###                    
### Tools: print, plot. summary.
###
### Copyright: Alejandro Jara, 2007
### Last modification: 05-07-2007.
###


"print.TDPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"summary.TDPdensity"<-function(object, hpd=TRUE, ...) 
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

### Parameters

    if(is.null(object$prior$a0))
    {
      dimen<-2
      coef.p<-object$coefficients[1:dimen]
      mat<-matrix(thetasave[,1:dimen],ncol=1)
    }
    else
    {
      dimen<-3
      coef.p<-object$coefficients
      mat<-thetasave[,1:dimen]

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

    ans$parms<-coef.table
    ans$nrec<-object$nrec

    class(ans) <- "summaryTDPdensity"
    return(ans)
}


"print.summaryTDPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$parms)) {
        cat("\nParameters:\n")
        print.default(format(x$parms, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}



"plot.TDPdensity"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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


   if(is(x, "TDPdensity"))
   {

      if(output=="density")
      {

      # Density estimation
	
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))

        title1<-paste("Density of",x$varnames[1],sep=' ')
	aa<-hist(x$y,plot=F,)
 	maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
	miny<-min(x$y)
	maxy<-max(x$y)
	deltay<-(maxy-miny)*0.2
	miny<-miny-deltay
	maxy<-maxy+deltay
	      
	hist(x$y,probability=T,xlim=c(min(x$grid),max(x$grid)),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")
        lines(x$grid,x$fun,lwd=2)
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
               if(pnames[i]=="ncluster"||pnames[i]=="k")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                 fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }
           
           if(is.null(x$prior$aa0))
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
            if(pnames[poss]=="ncluster"||pnames[poss]=="k")
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



