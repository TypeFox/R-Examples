### DPMdencens.R                    
### Fit a DPM of log-normal models for multivariate interval
### censored data.
###
### Copyright: Alejandro Jara, 2010-2012.
###
### Last modification: 09-07-2010.
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

"DPMdencens"<-
function(left,right,ngrid=100,grid=NULL,prior,mcmc,state,status)
UseMethod("DPMdencens")

"DPMdencens.default"<-
function(left,
         right, 
         ngrid=100,
         grid=NULL,
         prior,
         mcmc,
         state,
         status)
{
       #########################################################################################
       # data structure
       #########################################################################################

  	     cl <- match.call()

  	     nrec <- nrow(left)
  	     nvar <- ncol(left)
         if(is.null(nvar))
         { 
            nvar <- 1
            left <- as.matrix(left)
            right <- as.matrix(right)
            nrec <- nrow(left)
		 }

         y <- matrix(0,nrow=nrec,ncol=nvar)
         llower <- matrix(0,nrow=nrec,ncol=nvar)
         lupper <- matrix(0,nrow=nrec,ncol=nvar)
         tint <- matrix(0,nrow=nrec,ncol=nvar)

         for(i in 1:nvar)
         {
             tint[(!is.na(left[,i]) & is.na(right[,i])) ,i] <- 3
             tint[(!is.na(left[,i]) & !is.na(right[,i])) ,i] <- 2
			 tint[(is.na(left[,i]) & !is.na(right[,i])) ,i] <- 1
             tint[(is.na(left[,i]) & is.na(right[,i])) ,i] <- 4
             tint[(!is.na(left[,i]) & (left[,i]==right[,i])) ,i] <- 5
         }

         for(i in 1:nvar)
         {
			 llower[tint[,i]==2,i] <- left[tint[,i]==2,i]
             lupper[tint[,i]==2,i] <- right[tint[,i]==2,i]
             y[tint[,i]==2,i] <- 0.5*(left[tint[,i]==2,i]+right[tint[,i]==2,i])	

             lupper[tint[,i]==1,i] <- right[tint[,i]==1,i]
             y[tint[,i]==1,i] <- -1+right[tint[,i]==1,i]

             llower[tint[,i]==3,i] <- left[tint[,i]==3,i]
             y[tint[,i]==3,i] <- left[tint[,i]==3,i]+1
             
             y[tint[,i]==4,i] <- rnorm(1)
             
             y[tint[,i]==5,i] <- left[tint[,i]==5,i]
         }

       #########################################################################################
       # grid. Note: dimension of grid muxt be ngrid times nvar
       #########################################################################################

         if(is.null(grid))
         {
            grid <- matrix(0,nrow=ngrid,ncol=nvar)              
            for(i in 1:nvar)
            {
                tmp <- na.omit(c(left[,i],right[,i]))
                ll <- min(tmp)-1.0*sqrt(var(tmp))
				mm <- max(tmp)+1.0*sqrt(var(tmp))
                grid[,i] <- seq(ll,mm,len=ngrid)
            }
         }
		 else
         {
            grid <- as.matrix(grid)
            ngrid <- length(grid[,1])
         }

       #########################################################################################
       # change working directory (if requested..)
       #########################################################################################
       #  if(!is.null(work.dir))
       #  {
       #     cat("\n Changing working directory to ",work.dir,"\n")
       #     old.dir <- getwd()  # by default work in current working directory
       #     setwd(work.dir)
       #  }

	   #########################################################################################
	   # Prior specification
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
			psiinv1 <- matrix(var(y),nvar,nvar)
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
				m1[i] <- mean(y[,i])+rnorm(1,0,100)
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
         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave <- mcmc$nsave

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
				counter<-0
				for(j in 1:nvar)
				{
					muclus[i,j] <- m1[j]
					for(k in j:nvar)
					{
						counter <- counter+1
						sigmaclus[i,counter] <- psiinv1[j,k]
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
		 }    

	   #########################################################################################
	   # output
	   #########################################################################################
		 thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
         randsave <- matrix(0,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
         funi <- matrix(0,nrow=ngrid,ncol=nvar)
        
         nupp <- nvar*(nvar+1)/2 - nvar

	     ngridb <- nupp*ngrid*ngrid
         if(nvar==1) ngridb <- 1

         fbiv <- rep(0,ngridb)

	   #########################################################################################
	   # working space
	   #########################################################################################

         seed <- c(sample(1:29000,1),sample(1:29000,1))

		 muwork <- rep(0,nvar)
		 sigmawork <- matrix(0,nrow=nvar,ncol=nvar)
		 workm1 <- matrix(0,nrow=nvar,ncol=nvar)
		 workm2 <- matrix(0,nrow=nvar,ncol=nvar)
		 workm3 <- matrix(0,nrow=nvar,ncol=nvar)
		 workv1 <- rep(0,nvar)
		 workv2 <- rep(0,nvar)

		 ccluster <- rep(0,nrec)
		 cstrt <- matrix(0,nrow=nrec,ncol=nrec)
		 iflag <- rep(0,nvar)
		 prob <- rep(0,nrec+100)
		 ywork <- rep(0,nvar)
		 workmh1 <- rep(0,(nvar*(nvar+1)/2))
		 workmh2 <- rep(0,(nvar*(nvar+1)/2))

	   #########################################################################################
	   # calling the fortran code
	   #########################################################################################
         
         foo <- .Fortran("dpmdenscens",
					nrec      =as.integer(nrec),
                    nvar      =as.integer(nvar),
                    tint      =as.integer(tint),
                    llower    =as.double(llower),
                    lupper    =as.double(lupper),
                    ngrid     =as.integer(ngrid),
                    grid      =as.double(grid),
                    ngridb    =as.integer(ngridb),
					a0b0      =as.double(a0b0),
                    m1rand    =as.integer(m1rand),
                    nuvec     =as.integer(nuvec),
					psiinv2   =as.double(psiinv2),
					tau       =as.double(tau),
					s2inv     =as.double(s2inv),
					s2invm2   =as.double(s2invm2),
					mcmc      =as.integer(mcmcvec),
					nsave     =as.integer(nsave),
					alpha     =as.double(alpha),		
					k0		  =as.double(k0),
					m1        =as.double(m1),		
					muclus    =as.double(muclus),		 		
					ncluster  =as.integer(ncluster),
					psi1      =as.double(psi1),
					psiinv1   =as.double(psiinv1),
					ss        =as.integer(ss),
					sigmaclus =as.double(sigmaclus),
					y		  =as.double(y),
					randsave  =as.double(randsave),
					thetasave =as.double(thetasave),
                    fbiv      =as.double(fbiv),
                    funi      =as.double(funi),
					seed      =as.integer(seed),
					muwork    =as.double(muwork),
					sigmawork =as.double(sigmawork),
					workm1    =as.double(workm1),
					workm2    =as.double(workm2),
					workm3    =as.double(workm3),
					workv1    =as.double(workv1),
					workv2    =as.double(workv2),
					ccluster  =as.integer(ccluster),
					cstrt     =as.integer(cstrt),
					iflag     =as.integer(iflag),
					prob      =as.double(prob),
					workmh1   =as.double(workmh1),
					workmh2   =as.double(workmh2),
					ywork     =as.double(ywork),
					PACKAGE="DPpackage")	

	   #########################################################################################
	   # save state
	   #########################################################################################

         model.name <- "DPM model for interval-censored data"

         varnames <- paste("var",1:nvar,sep="")

         state <- list(
					   alpha=foo$alpha,
					   m1=matrix(foo$m1,nrow=nvar,ncol=1),
                       muclus=matrix(foo$muclus,nrow=nrec+100,ncol=nvar),
					   ncluster=foo$ncluster,
                       psi1=matrix(foo$psi1,nrow=nvar,ncol=nvar),
                       k0=foo$k0,
                       sigmaclus=matrix(foo$sigmaclus,nrow=nrec+100,ncol=nuniq),
                       ss=foo$ss,
                       y=matrix(foo$y,nrow=nrec,ncol=nvar)
                       )

		 randsave <- matrix(foo$randsave,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)


         pnames <- paste("m1",varnames,sep="-")
         pnames <- c(pnames,"k0")
         for(i in 1:nvar)
         {  
             for(j in i:nvar)
             {
                 tmp1 <- paste("var",i,sep="") 
                 tmp2 <- paste("var",j,sep="")
                 tmp3 <- paste(tmp1,tmp2,sep=":")
                 tmp4 <- paste("psi1",tmp3,sep="-")
                 pnames <- c(pnames,tmp4)
			 }
         }
         pnames <- c(pnames,"ncluster","alpha")
         colnames(thetasave) <- pnames
 
         coeff <- apply(thetasave,2,mean)

         
         tmp1 <- paste("mu",1:nvar,sep="")   
         tmp2 <- paste("sigma",1:(nvar*(nvar+1)/2),sep="")   
         tmp3 <- paste("id=",1:nrec,sep="")

         pnamesre <- paste(rep(c(tmp1,tmp2),nrec),rep(paste("id",1:nrec,sep="="),each=(nvar+nvar*(nvar+1)/2)),sep=".")
         pnamesre <- c(pnamesre,paste(c(tmp1,tmp2),"pred",sep="."),paste(paste("var",1:nvar,sep=""),"pred",sep="."))

         colnames(randsave) <- pnamesre

		 save.state <- list(thetasave=thetasave,randsave=randsave)

         ff <- matrix(foo$funi,nrow=ngrid,ncol=nvar)

         funi <- NULL
         for(i in 1:nvar)
         {
             funi[[i]] <- ff[,i]
         }

		 fbiv <- NULL
         if(nvar > 1)
         {
            count <- 0
            beg <- 1
            end <- ngrid*ngrid 
            for(i in 1:(nvar-1))
            {
                for(j in (i+1):nvar)
                {
                   count <- count+1
                   fbiv[[count]] <- matrix(foo$fbiv[beg:end],nrow=ngrid,ncol=ngrid, byrow = TRUE)
                   beg <- end+1
                   end <- end + ngrid*ngrid
                }
            }  
         }

		 z <- list(modelname=model.name,
				   call=cl,
                   coefficients=coeff,
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   save.state=save.state,
                   nrec=nrec,
                   nvar=nvar,
                   funi=funi,
                   fbiv=fbiv, 
                   grid=grid,
                   varnames=varnames)

	 cat("\n\n")
	 class(z)<-c("DPMdencens")
	 z 
}

###
### Tools for DPMdencens: print, summary, plot
###
### Copyright: Alejandro Jara, 2010.
### Last modification: 19-07-2010.


"print.DPMdencens"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)

    cat("\nPosterior Inference of Parameters:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Variables:",x$nvar,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.DPMdencens"<-function(object, hpd=TRUE, ...) 
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

### Baseline Information

    mat <- NULL
    coef.p <- NULL
    
    dimen1 <- object$nvar+object$nvar*(object$nvar+1)/2+1
    
    for(i in 1:dimen1)
    {
	    coef.p<-c(coef.p,object$coefficients[i])
		mat <- cbind(mat,thetasave[,i])
    }

    if(dimen1>0)
    {
    
       coef.m <-apply(mat, 2, median)    
       coef.sd<-apply(mat, 2, sd)
       coef.se<-apply(mat, 2, stde)

       if(hpd)
       {             
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
    }

### Precision parameter

    dimen1 <- object$nvar+object$nvar*(object$nvar+1)/2+1
    
    if(is.null(object$prior$a0))
    {
      dimen2 <- 1
      coef.p<-object$coefficients[(dimen1+1)]
      mat<-matrix(thetasave[,(dimen1+1)],ncol=1)
    }
    else
    {
      dimen2 <- 2
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

    class(ans) <- "summaryDPMdencens"
    return(ans)
}


"print.summaryDPMdencens"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    	     
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




"plot.DPMdencens"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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


   if(is(x, "DPMdencens"))
   {

      if(output=="density")
      {

      # Density estimation
	
    	par(ask = ask)
	    layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))
		start<-(x$nrec+1)*x$nvar+(x$nrec+1)*x$nvar*(x$nvar+1)/2

        if(x$nvar==1)
        {
           title1 <- paste("Density of",x$varnames[1],sep=' ')
           plot(x$grid[,1],x$funi[[1]],type="l",lwd=2,main=title1,xlab="values", ylab="density")
        }
        else
        {
           for(i in 1:x$nvar)
           {
               title1 <- paste("Density of",x$varnames[i],sep=' ')
			   plot(x$grid[,i],x$funi[[i]],type="l",lwd=2,main=title1,xlab="values", ylab="density")
           }
           
           count <- 0
           for(i in 1:(x$nvar-1))
           {
               for(j in (i+1):x$nvar)
               {
                   count <- count +1
                   varsn <- paste(x$varnames[i],x$varnames[j],sep="-")
                   title1 <- paste("Density of ",varsn,sep='')
                   xx <- matrix(x$grid[,i],ncol=1)
                   yy <- matrix(x$grid[,j],ncol=1)
                   z <- x$fbiv[[count]] 
                   colnames(xx) <- x$varnames[i]
                   colnames(yy) <- x$varnames[j]
 
				   contour(xx,yy,z,main=title1,xlab=x$varnames[i],ylab=x$varnames[j])
                   persp(xx,yy,z,xlab=x$varnames[i],ylab=x$varnames[j],zlab="density",theta=-30,phi=15,expand = 0.9, ltheta = 120,main=title1)
                }
		   }
        }
      }
      else
      {

        if(is.null(param))
        {
           pnames <- colnames(x$save.state$thetasave)
           n <- dim(x$save.state$thetasave)[2]
           cnames <- names(x$coefficients)
           
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
            pnames <- colnames(x$save.state$thetasave)
            n <- dim(x$save.state$thetasave)[2]
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






