### LDPDdoublyint.R                    
### Fit a semiparametric aft model for multivariate doubly interval
### censored survival data.
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 23-01-2010.
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

"LDPDdoublyint"<-
function(onset,failure,p,xonset,q,xfailure,xpred,grid,prior,mcmc,state,status,work.dir=NULL)
UseMethod("LDPDdoublyint")

"LDPDdoublyint.default"<-
function(onset,
         failure, 
         p,
         xonset,
         q,
         xfailure,
         xpred,
         grid,
         prior,
         mcmc,
         state,
         status,
         work.dir=NULL)
{

       #########################################################################################
       # data structure
       #########################################################################################

  	     cl <- match.call()

  	     nsubject <- nrow(onset)
  	     nvar <- ncol(onset)/2

         x <- cbind(xonset,xfailure)  	 
         
         y <- NULL
         llower <- NULL
         lupper <- NULL
         j <- 0
         for(i in 1:nvar)
         {
            j <- j+1
            llower <- cbind(llower,onset[,j])
            y <- cbind(y,log((onset[,j]+onset[,j+1])/2))
            j <- j+1
            lupper <- cbind(lupper,onset[,j])
         }

         j <- 0
         for(i in 1:nvar)
         {
            j <- j+1         
            llower <- cbind(llower,failure[,j])
            tmp <- ((failure[,j]+failure[,j+1])/2)-((onset[,j]+onset[,j+1])/2)
            tmp[tmp==0] <- 0.01
            y <- cbind(y,log(tmp))
            j <- j+1
            lupper <- cbind(lupper,failure[,j])
            
         }
 
         npred <- nrow(xpred)
         ngrid <- ncol(grid)

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
	   # Prior specification
	   #########################################################################################

         maxm <- prior$maxm
         
         nu <- c(prior$nu,prior$nub)

         if(is.null(prior$a0))
         {
            a0b0 <- c(-1,-1,0)
            ap <- prior$a
            apr <- 0
         }
         else
         {
            a0b0 <- c(prior$a0,prior$b0,prior$q)
            ap <- 0
            apr <- 1
            if(prior$a0<0 || prior$b0<0)
            { 
                   stop("The parameters of the Beta prior for the a parameter must be possitive.\n")     
            }

            if(prior$q<0 || prior$q>1)
            { 
                   stop("The mixture proportion must be in (0,1).\n")     
            }
         }

         if(is.null(prior$mub))
         {
            a0b0 <- c(a0b0,0,-1)
            bp <- prior$b
            bpr <- 0
         }
         else
         {
            a0b0 <- c(a0b0,prior$mub,prior$sigmab)
            bpr <- 1
            bp <- abs(rnorm(1,mean=prior$mub,sd=sqrt(prior$sigmab)))
         }
         
         alpha <- c(ap,bp)
           
         tinv1 <- prior$tinv
         tinv2 <- prior$tbinv
          
         psiinv <- solve(prior$S0)
         smu <- psiinv%*%prior$m0
           
	   #########################################################################################
	   # mcmc specification
	   #########################################################################################

         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave <- mcmc$nsave

         if(is.null(mcmc$tune1))
         {  
            tune1 <- 1.1
         }
         else
         {
            tune1 <- mcmc$tune1
         }

         if(is.null(mcmc$tune2))
         {
            tune2 <- 1.1
         }
         else
         {
            tune2 <- mcmc$tune2
         }
         a0b0 <- c(a0b0,tune1,tune2)

	   #########################################################################################
	   # Starting values
	   #########################################################################################

         sigma <- matrix(0,nrow=2*nvar,ncol=2*nvar)

         beta <- NULL 
         sigma <- NULL

         for(i in 1:nvar)
         {
             startp <- (i-1)*p+1
             endp<- (i-1)*p+p
             
             fit0 <- glm.fit(xonset[,(startp:endp)], y[,i], family=gaussian())   
             betatmp <- coefficients(fit0)
             beta <- c(beta,betatmp)
             s2 <- var(residuals(fit0))
             sigma <- c(sigma,s2)
         }
         for(i in 1:nvar)
         {
             startp <- (i-1)*q+1
             endp <- (i-1)*q+q
             fit0 <- glm.fit(xfailure[,(startp:endp)], y[,nvar+i], family=gaussian())   
             betatmp <- coefficients(fit0)
             beta <- c(beta,betatmp)
             s2 <- var(residuals(fit0))
             sigma <- c(sigma,s2)
         }

         sigma <- diag(sigma,2*nvar)
         mub <- beta
         sigmab <- diag(1,nvar*(p+q))

	   #########################################################################################
	   # output
	   #########################################################################################
         acrate <- rep(0,2)
         npar <- (2*nvar*(2*nvar+1)/2)+nvar*(p+q)+(nvar*(p+q)*(nvar*(p+q)+1)/2)+3

         f <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)
         h <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)
         thetasave <- matrix(0,nrow=nsave,ncol=npar)
         randsave <- matrix(0,nrow=nsave,ncol=npred*2*nvar)
         
	   #########################################################################################
	   # parameters depending on status
	   #########################################################################################

	     if(status==TRUE)
	     {
		    sigmainv <- solve(sigma) 
		    sigmabinv <- solve(sigmab)
 	        b <- matrix(0,nrow=maxm,ncol=nvar*(p+q))
 	        for(i in 1:1)
 	        {
 	           b[i,] <- beta
 	        }
 	        for(i in 2:maxm)
 	        {
 	           b[i,] <- rep(0,nvar*(p+q))
 	        }

 	        ss <- rep(1,nsubject)
 	        ncluster <- 1
	     }	
	     if(status==FALSE)
	     {
	        alpha <- state$alpha 
		    b <- state$b
		    sigma <- state$sigma
	     	sigmainv <- solve(state$sigma)
		    mub <- state$mb
		    sigmab <- state$Sb
		    sigmabinv <- solve(state$Sb)
			y <- state$y
		    ss <- state$ss
		    ncluster <- state$ncluster
	     }

	   #########################################################################################
	   # working space
	   #########################################################################################

         cstrt <- matrix(0,nrow=maxm,ncol=nsubject)
         ccluster <- rep(0,maxm)
         prob <- rep(0,maxm)
         prob2 <- rep(0,maxm)
         weights <- rep(0,maxm)

         iflagc <- rep(0,nvar*(p+q))
         theta <- rep(0,nvar*(p+q))
         workmc <- matrix(0,nrow=nvar*(p+q),ncol=nvar*(p+q))
         workmc2 <- matrix(0,nrow=nvar*(p+q),ncol=nvar*(p+q))
         workmc3 <- matrix(0,nrow=nvar*(p+q),ncol=nvar*(p+q))
         workmhc <- rep(0,nvar*(p+q)*(nvar*(p+q)+1)/2)
         workmhc2 <- rep(0,nvar*(p+q)*(nvar*(p+q)+1)/2)
         workvc <- rep(0,nvar*(p+q))         

         iflagn <- rep(0,2*nvar)
         workmn <- matrix(0,nrow=2*nvar,ncol=2*nvar)
         workmn2 <- matrix(0,nrow=2*nvar,ncol=2*nvar)
         workmn3 <- matrix(0,nrow=2*nvar,ncol=2*nvar)
         workmhn <- rep(0,2*nvar*(2*nvar+1)/2)
         workmhn2 <- rep(0,2*nvar*(2*nvar+1)/2)
         workvn <- rep(0,2*nvar)
         workvn2 <- rep(0,2*nvar)
         workvn3 <- rep(0,2*nvar)
         ztz <- matrix(0,nrow=nvar*(p+q),ncol=nvar*(p+q))
         zty <- rep(0,nvar*(p+q))

	     seed <- c(sample(1:29000,1),sample(1:29000,1))

         model <- matrix(0,nrow=2*nvar,ncol=(p+q)) 
         possi <- matrix(0,nrow=2*nvar,ncol=(p+q)) 
         
         fw <- rep(0,ngrid)
         fw2 <- rep(0,ngrid)

         varind <- rep(1,p)
         for(i in 2:nvar)
         {
             varind <- c(varind,rep(i,p))
         }
         for(i in 1:nvar)
         {
             varind <- c(varind,rep(nvar+i,q))
		 }

	   #########################################################################################
	   # calling the fortran code
	   #########################################################################################
         
         foo <- .Fortran("ldpddoublyintsba",
					nsubject  =as.integer(nsubject),
					nvar      =as.integer(nvar),
					p         =as.integer(p),
					q         =as.integer(q),
					ngrid     =as.integer(ngrid),
					npred     =as.integer(npred),	 	
					x         =as.double(x),
					llower    =as.double(llower),
					lupper    =as.double(lupper),
					grid      =as.double(grid),
					xpred     =as.double(xpred),
					maxm      =as.integer(maxm),
					a0b0      =as.double(a0b0),
					nu        =as.integer(nu),
					tinv1     =as.double(tinv1),
					smu       =as.double(smu),
					psiinv    =as.double(psiinv),
					tinv2     =as.double(tinv2),
					mcmc      =as.integer(mcmcvec),
					nsave     =as.integer(nsave),
					ncluster  =as.integer(ncluster),
					ss        =as.integer(ss),
					alpha     =as.double(alpha),
					b         =as.double(b),
					sigma     =as.double(sigma),
					sigmainv  =as.double(sigmainv),
					mub       =as.double(mub),
					sigmab    =as.double(sigmab),
					sigmabinv =as.double(sigmabinv),
					y         =as.double(y),
					acrate    =as.double(acrate),
					f         =as.double(f),
					h         =as.double(h),
					thetasave =as.double(thetasave),
					randsave  =as.double(randsave),
					seed      =as.integer(seed),
					model     =as.integer(model),
					possi     =as.integer(possi),
                    varind    =as.integer(varind),
					cstrt     =as.integer(cstrt),
					ccluster  =as.integer(ccluster),
					prob      =as.double(prob),
					prob2     =as.double(prob2),
					weights   =as.double(weights),
					iflagc    =as.integer(iflagc),
					theta     =as.double(theta),
					workmc    =as.double(workmc),
					workmc2   =as.double(workmc2),
					workmc3   =as.double(workmc3),
					workmhc   =as.double(workmhc),
					workmhc2  =as.double(workmhc2),
					workvc    =as.double(workvc),
					iflagn    =as.integer(iflagn),
					workmn    =as.double(workmn),
					workmn2   =as.double(workmn2),
					workmn3   =as.double(workmn3),
					workmhn   =as.double(workmhn),
					workmhn2  =as.double(workmhn2),
					workvn    =as.double(workvn),
					workvn2   =as.double(workvn2),
					workvn3   =as.double(workvn3),
					ztz       =as.double(ztz),
					zty       =as.double(zty),
					fw        =as.double(fw),
					fw2       =as.double(fw2),
					PACKAGE="DPpackage")	

	   #########################################################################################
	   # save state
	   #########################################################################################

         model.name<-"Linear Dependent Poisson-Dirichlet Survival Model for Doubly-Interval Censored Data"

         nametmp<-c(paste("onset",1:nvar,sep=""),
                    paste("time",1:nvar,sep="")
                    ) 
 
         pnames<-NULL
         for(i in 1:(2*nvar))
         {
            for(j in i:(2*nvar))
            {
                tmp<-paste(nametmp[i],nametmp[j],sep=";")
                tmp<-paste("Sigma (",tmp,")",sep="")
                pnames<-c(pnames,tmp)
            }
         }

         nametmp <- colnames(x) 
         pnames <- c(pnames,paste("mb",nametmp,sep="-"))

         for(i in 1:(nvar*(p+q)))
         {
            for(j in i:(nvar*(p+q)))
            {
                tmp<-paste(nametmp[i],nametmp[j],sep=";")
                tmp<-paste("Sb (",tmp,")",sep="")
                pnames<-c(pnames,tmp)
            }
         }
         
         pnames <- c(pnames,"ncluster","a","b")
         
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=npar)
         colnames(thetasave) <- pnames
         
         coeff <- apply(thetasave,2,mean)		


         nametmp<-c(paste("onset",1:nvar,sep=""),
                    paste("time",1:nvar,sep="")
                    ) 
         nametmp <- rep(nametmp,npred)           
         randsave <- matrix(foo$randsave,nrow=nsave,ncol=npred*2*nvar)
         colnames(randsave)<-nametmp
         
         coeffr <- apply(randsave,2,mean)		
         
         medtable <- matrix(coeffr,nrow=npred,ncol=2*nvar,byrow=T)
         nametmp<-c(paste("onset",1:nvar,sep=""),
                    paste("time",1:nvar,sep="")
                    ) 
         colnames(medtable)<-nametmp
         nametmp<-paste("prediction",1:npred)
         rownames(medtable)<-nametmp

		 save.state <- list(thetasave=thetasave,
							randsave=randsave)
	                    

		 state <- list(	alpha=foo$alpha, 
						b=matrix(foo$b,nrow=maxm,ncol=nvar*(p+q)),
						sigma=matrix(foo$sigma,nrow=2*nvar,ncol=2*nvar),
						mb=foo$mub,
						Sb=matrix(foo$sigmab,nrow=nvar*(p+q),ncol=nvar*(p+q)),
						y=matrix(foo$y,nrow=nsubject,ncol=2*nvar),
						ss=foo$ss,
						ncluster=foo$ncluster)
 
         acrate <- foo$acrate[1]
         dptest <- foo$acrate[2]


		 z<-list(modelname=model.name,
				 coefficients=coeff,
				 call=cl,
				 acrate=acrate,
				 dptest=dptest,
                 prior=prior,
                 mcmc=mcmc,
                 state=state,
                 save.state=save.state,
                 nsubject=foo$nsubject,
                 nsave=nsave,
                 p=foo$p,
                 q=foo$q,
                 nvar=2*foo$nvar,
                 npred=foo$npred,
                 ngrid=ngrid,
                 grid=grid,
                 f=matrix(foo$f,nrow=npred*2*nvar,ncol=ngrid),
                 h=matrix(foo$h,nrow=npred*2*nvar,ncol=ngrid),
                 medtable=medtable,
                 acrate=foo$acrate,
                 work.dir=work.dir)

	 cat("\n\n")
	 class(z)<-c("LDPDdoublyint")
	 z 
}

###
### Tools for LDPDdoublyint: print, summary, plot, predict
###
### Copyright: Alejandro Jara, 2007
### Last modification: 08-13-2007.



"print.LDPDdoublyint"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)

    cat("\nP(DP|data)=",x$dptest,"\n")

	cat("\nPosterior Inference of Parameters:\n")
	print.default(format(x$coefficients, digits = digits), print.gap = 2, 
	quote = FALSE)

    if(!is.null(x$acrate)) 
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate)    
    }   
    cat("\n\nNumber of Subjects:",x$nsubject)
    cat("\nNumber of Variables:",x$nvar)
    cat("\n\n")
    invisible(x)
}



"summary.LDPDdoublyint"<-function(object, hpd=TRUE, ...) 
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

### Kernel 

	nvar <- object$nvar

    dimen1 <- nvar*(nvar+1)/2
    
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

       ans$kernel <- coef.table
    }


### Baseline

	p <- object$p
	q <- object$q

    nt <- (p+q)*(nvar/2)

    dimen2 <- nt + nt*(nt+1)/2
    
    coef.p <- object$coefficients[(dimen1+1):(dimen1+dimen2)]
    mat <- thetasave[,(dimen1+1):(dimen1+dimen2)]

    if(dimen2>0)
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

       ans$base <- coef.table
    }


### Precision 

    dimen3 <- 3

    coef.p <- object$coefficients[(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]
    mat <- thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]

    if(dimen3>0)
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

       ans$prec <- coef.table
    }

    ans$nsubject <- object$nsubject
    ans$nvar <- object$nvar
    ans$acrate <- object$acrate
    ans$dptest <- object$dptest

    class(ans) <- "summaryLDPDdoublyint"
    return(ans)
}


"print.summaryLDPDdoublyint"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    	     
    cat("\nP(DP|data)=",x$dptest,"\n")

	cat("\nKernel covariance matrix:\n")
	print.default(format(x$kernel, digits = digits), print.gap = 2, 
	quote = FALSE)

	cat("\nCentering parameters:\n")
	print.default(format(x$base, digits = digits), print.gap = 2, 
	quote = FALSE)

	cat("\nPrecision parameters:\n")
	print.default(format(x$prec, digits = digits), print.gap = 2, 
	quote = FALSE)

    if(!is.null(x$acrate)) 
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate)    
    }   
    cat("\n\nNumber of Subject:",x$nsubject)
    cat("\nNumber of Variable:",x$nvar)
    cat("\n\n")
    invisible(x)
}


"plot.LDPDdoublyint"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
{

fancydensplot<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
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


    if(is(x, "LDPDdoublyint"))
    {
        if(is.null(param))
	    {
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:(n-2)){
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster" || pnames[i]=="k")
	       {
	          hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
	       }
	       else
	       {
                  fancydensplot(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }   
           }
           
           if(!is.null(x$prior$a0))
           {
               title1<-paste("Trace of",pnames[n-1],sep=" ")
               title2<-paste("Density of",pnames[n-1],sep=" ")       
               plot(x$save.state$thetasave[,n-1],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,n-1],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           if(!is.null(x$prior$mub))
           {
               title1<-paste("Trace of",pnames[n],sep=" ")
               title2<-paste("Density of",pnames[n],sep=" ")       
               plot(x$save.state$thetasave[,n],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
        }
        else
        {
            coef.p <- x$coefficients
			n <- length(coef.p)
			pnames<-names(coef.p)
			poss<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }

            if (poss==0 && param !="predictive") 
			{
				stop("This parameter is not present in the original model.\n")
			}

            if (param !="ncluster") 
            {
				par(ask = ask)
				layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
				title1<-paste("Trace of",pnames[poss],sep=" ")
				title2<-paste("Density of",pnames[poss],sep=" ")       
				plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
				fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            else
            {
				par(ask = ask)
				layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
				title1<-paste("Trace of",pnames[poss],sep=" ")
				title2<-paste("Density of",pnames[poss],sep=" ")       
				plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
				hist(x$save.state$thetasave[,poss],main=title2,xlab="values", ylab="probability",probability=TRUE)
            }
        }
   }
}



"predict.LDPDdoublyint"<-
function(object,grid,compute.band=FALSE,hpd=TRUE, ...)
{
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


   if(is(object, "LDPDdoublyint"))
   {

	  ngrid <- object$ngrid
	  grid <- object$grid
	  nvar <- object$nvar/2
	  npred <- object$npred

	  f <- object$f
	  h <- object$h
	  med <- object$medtable

      if(compute.band)
	  {
		nsave <- object$nsave
		alpha <- 0.05
		tint <- 0
		if(hpd)tint <- 1
		llower <- NULL
		lupper <- NULL
      
		llower <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)
		lupper <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)

		llower2 <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)
		lupper2 <- matrix(0,nrow=npred*2*nvar,ncol=ngrid)

		workv1 <- rep(0,nsave)
		workv2 <- rep(0,ngrid)

		foo <- .Fortran("hpdspy",
            	nsave     =as.integer(nsave),
            	npred     =as.integer(npred),
				nvar      =as.integer(nvar),
                ngrid     =as.integer(ngrid),
                alpha     =as.double(alpha),
                tint      =as.integer(tint),
                workv1    =as.double(workv1), 
                workv2    =as.double(workv2),
                llower    =as.double(llower),
                lupper    =as.double(lupper),
                llower2   =as.double(llower2),
                lupper2   =as.double(lupper2),
				PACKAGE="DPpackage")	

		llower <- matrix(foo$llower,nrow=npred*2*nvar,ncol=ngrid)
		lupper <- matrix(foo$lupper,nrow=npred*2*nvar,ncol=ngrid)
		llower2 <- matrix(foo$llower2,nrow=npred*2*nvar,ncol=ngrid)
		lupper2 <- matrix(foo$lupper2,nrow=npred*2*nvar,ncol=ngrid)

        if(hpd)
		{
			limm <- apply(object$save.state$randsave, 2, hpdf)
			coef.l <- limm[1,]
			coef.u <- limm[2,]
        }
        else
        {
			limm <- apply(object$save.state$randsave, 2, pdf)
			coef.l <- limm[1,]
			coef.u <- limm[2,]
		}

		med.l <- matrix(coef.l,nrow=npred,ncol=2*nvar,byrow=T)
		med.u <- matrix(coef.u,nrow=npred,ncol=2*nvar,byrow=T)

		colnames(med.l) <- colnames(med)
		colnames(med.u) <- colnames(med)
		rownames(med.l) <- rownames(med)
		rownames(med.u) <- rownames(med)

      }
	  else
	  {
		llower <- NULL
		lupper <- NULL
		llower2 <- NULL
		lupper2 <- NULL
		med.l <- NULL
		med.u <- NULL
	  }


      out <- list(f=f,
				  h=h,
				  med=med,
				  med.l=med.l,
				  med.u=med.u,
				  grid=grid,
				  ngrid=ngrid,
				  nvar=nvar,
				  npred=npred,
				  llower=llower,
				  lupper=lupper,
				  llower2=llower2,
				  lupper2=lupper2,
				  modelname=object$modelname,
				  call=object$call,
				  compute.band=compute.band)


      class(out) <- c("predict.LDPDdoublyint")
      out
   }
}   


"print.predict.LDPDdoublyint"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)

	cat("\nPosterior Mean of Median Survival for Predictions:\n")
	print.default(format(x$med, digits = digits), print.gap = 2, 
	quote = FALSE)

	if(!is.null(x$med.l)) 
	{
		cat("\nLower CI Limit for the Median Survival for Predictions:\n")
		print.default(format(x$med.l, digits = digits), print.gap = 2, 
		quote = FALSE)
	}   

	if(!is.null(x$med.u)) 
	{
		cat("\nUpper CI Limit for the Median Survival for Predictions:\n")
		print.default(format(x$med.u, digits = digits), print.gap = 2, 
		quote = FALSE)
	}   

    cat("\n\n")
    invisible(x)
}


plot.predict.LDPDdoublyint<-function(x,ask=TRUE,xlim=NULL,nfigr=1,nfigc=1, ...)
{
    if(is(x, "predict.LDPDdoublyint"))
    {	
		par(ask = ask)
		layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
		compute.band <- x$compute.band

		if(compute.band)
		{

			k <- 0
			for(i in 1:x$npred)
			{
				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Onset",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1

					x1 <- x$grid[j,]
					y1 <- x$llower[k,]
					x2 <- x$grid[j,]
					y2 <- x$lupper[k,]
					aa <- rbind(x2,y2)[, order(-x2, y2)]
					x2 <- aa[1,]
					y2 <- aa[2,]

					plot(x$grid[j,],x$f[k,],type="l",lty=1,lwd=3,ylim=c(0,1),
						 ylab="Survival",xlab="Time",main=tmp3,bty="n")           
					polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
					lines(x$grid[j,],x$f[k,],lwd=3,lty=1)
				}

				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Failure",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1

					x1 <- x$grid[j,]
					y1 <- x$llower[k,]
					x2 <- x$grid[j,]
					y2 <- x$lupper[k,]
					aa <- rbind(x2,y2)[, order(-x2, y2)]
					x2 <- aa[1,]
					y2 <- aa[2,]

					plot(x$grid[j,],x$f[k,],type="l",lty=1,lwd=3,ylim=c(0,1),
						 ylab="Survival",xlab="Time",main=tmp3,bty="n")           
					polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
					lines(x$grid[j,],x$f[k,],lwd=3,lty=1)
				}
			}


			k <- 0
			for(i in 1:x$npred)
			{
				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Onset",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1

					x1 <- grid[j,]
					y1 <- x$llower2[k,]
					x2 <- grid[j,]
					y2 <- x$lupper2[k,]
					aa <- rbind(x2,y2)[, order(-x2, y2)]
					x2 <- aa[1,]
					y2 <- aa[2,]

					plot(x$grid[j,],x$lupper2[k,],type="l",lty=1,lwd=2,
						ylab="Hazard",xlab="Time",main=tmp3,bty="n",col="lightgray") 
					polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
					lines(x$grid[j,],x$h[k,],lwd=3,lty=1)
				}

				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Failure",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1

					x1 <- grid[j,]
					y1 <- x$llower2[k,]
					x2 <- grid[j,]
					y2 <- x$lupper2[k,]
					aa <- rbind(x2,y2)[, order(-x2, y2)]
					x2 <- aa[1,]
					y2 <- aa[2,]

					plot(x$grid[j,],x$lupper2[k,],type="l",lty=1,lwd=2,
						ylab="Hazard",xlab="Time",main=tmp3,bty="n",col="lightgray") 
					polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
					lines(x$grid[j,],x$h[k,],lwd=3,lty=1)
				}
			}

		}
		else
		{

			k <- 0
			for(i in 1:x$npred)
			{
				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Onset",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1
					plot(x$grid[j,],x$f[k,],type="l",lty=1,lwd=3,ylim=c(0,1),
						ylab="Survival",xlab="Time",main=tmp3,bty="n")           
				}

				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Failure",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1
					plot(x$grid[j,],x$f[k,],type="l",lty=1,lwd=3,ylim=c(0,1),
							ylab="Survival",xlab="Time",main=tmp3,bty="n")           
				}
			}
   
			k <- 0
			for(i in 1:x$npred)
			{
				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Onset",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1
					plot(x$grid[j,],x$h[k,],type="l",lty=1,lwd=3,
						ylab="Hazard",xlab="Time",main=tmp3,bty="n")           
				}

				for(j in 1:x$nvar)
				{
					tmp1 <- paste("Prediction",i,sep=":")
					tmp2 <- paste("Failure",j,sep=":")
					tmp3 <- paste(tmp1,tmp2,sep="-")
					k <- k+1
					plot(x$grid[j,],x$h[k,],type="l",lty=1,lwd=3,
						ylab="Hazard",xlab="Time",main=tmp3,bty="n")           
				}
			}
		}
    }
}    

