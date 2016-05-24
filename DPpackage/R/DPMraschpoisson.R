### DPMraschpoisson.R                   
### Fit a Rasch Poisson model with a Dirichlet mixture of normals prior
### for the random effect distribution
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 04-09-2009.
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


DPMraschpoisson <- function(y,prior,mcmc,offset=NULL,state,status,
                            grid=seq(-10,10,length=1000),data=sys.frame(sys.parent()),compute.band=FALSE)
UseMethod("DPMraschpoisson")

DPMraschpoisson.default <-
function(y,
         prior,
         mcmc,
         offset=NULL,         
         state,
         status,
         grid=seq(-10,10,length=1000),
         data=sys.frame(sys.parent()),
         compute.band=FALSE)
{
         #########################################################################################
         # call parameters
         #########################################################################################
           cl <- match.call()
     	   y <- as.matrix(y)
	  
         #########################################################################################
         # data structure
         #########################################################################################
		   nsubject <- nrow(y)
	       p <- ncol(y)

           ywork <- y
          
           datastrm <- NULL
           nmissing <- 0
           total <- 0
          
           for(i in 1:nsubject)
           {
               for(j in 1:p)
               {
                  if(is.na(y[i,j]))
                  {
                     nmissing <- nmissing+1
                     datastrm <- rbind(datastrm,c(i,j))   
                  }
                  else
                  {
                     total<-total+y[i,j]            
                  }
                  
               }
		   }
          
           nrec <- nsubject*p-nmissing
          
           if(nmissing>0)
           {
              imiss <- 1 
              for(i in 1:nmissing)
              {
                   ywork[datastrm[i,1],datastrm[i,2]] <- rpois(1,total/nrec)               
              }
           }
           else
           {
              imiss <- 0
              nmissing <- 1
              datastrm <- matrix(0,nrow=1,ncol=2)
           }


         #########################################################################################
         # prior information
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

           tauk1 <- prior$tauk1
           if(tauk1 < 0)
           { 
			  stop("The parameter of the Gamma prior for the kernel variance must be possitive.\n")     
           }


           homo.var <- 0
           if(is.null(prior$tauk2))
           {
              taus1 <- prior$taus1
              taus2 <- prior$taus2
              if(taus1 < 0 || taus2 < 0)
              { 
                  stop("The parameters of the Gamma prior for the gamma centering distribution must be possitive.\n")     
              }
			  tauk2 <- 2.01
              tauk2rand <- 1
           } 
           else
           {
			  taus1 <- -1
              taus2 <- -1
              tauk2 <- prior$tauk2
              tauk2rand <- 0
           }

  	       if(is.null(prior$taub1))
  	       {
              taub1 <- -1
              taub2 <- -1
              sigmab <- prior$sigmab
              sigmarand <- 0
  	       }
  	       else
  	       {
              taub1 <- prior$taub1
              taub2 <- prior$taub2
              sigmab <- 1
  	          sigmarand <- 1
		    }

			tau <- c(tauk1,taub1,taub2,taus1,taus2)

			if(is.null(prior$m0))
  	        {
  	            s0 <- -1
  	            m0 <- 0
  	            mub <- prior$mub
  	            murand <- 0
  	        }
  	        else
  	        {
  	            s0 <- prior$s0
				m0 <- prior$m0
       	        mub <- rnorm(1,m0,sqrt(s0))
                murand <- 1
            }     


            if(is.null(prior$N))
            {
                maxn <- 50
            }
            else
            {
                maxn <- prior$N
			}


            b0 <- prior$beta0
            prec <- solve(prior$Sbeta0)
            sb <- cbind(prec%*%b0,b0)

            if(length(b0)!=(p-1))
            { 
			   stop("Error in the dimension of the mean of the normal prior for the difficulty parameters.\n")     
            }

            if(nrow(prec)!=(p-1) || ncol(prec)!=(p-1))
			{ 
			   stop("Error in the dimension of the covariance of the normal prior for the difficulty parameters.\n")     
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
              nsave <- mcmc$nsave
           }


         #########################################################################################
         # output
         #########################################################################################
		   acrate <- 0
		   cpo <- matrix(0,nrow=nsubject,ncol=p)
		   cpov <- rep(0,nsubject)

		   ngrid <- length(grid)

           thetasave <- matrix(0,nrow=nsave,ncol=p+6)
           randsave <- matrix(0,nrow=nsave,ncol=nsubject+1)
           densave <- matrix(0,nrow=nsave,ncol=ngrid)
           cdfsave <- matrix(0,nrow=nsave,ncol=ngrid)

         #########################################################################################
         # MLE estimation
         #########################################################################################
         
           RaschMLE <- function(y,nitem,nsubject,offset)
           {
             ywork2 <- NULL
             roffset <- NULL
             id <- NULL
             x <- NULL
             count <- 0
             for(i in 1:nsubject)
             {
				 ywork2 <- c(ywork2,y[i,])
				 roffset <- c(roffset,offset[i,])
				 id <- c(id,rep(i,nitem))

                 aa <- diag(-1,nitem)
				 aa[,1] <- 1
                 x <- rbind(x,aa)
             }
             out <- NULL
             library(nlme)
			 library(MASS)
             fit0 <- glmmPQL(ywork2~x-1+offset(roffset),random = ~ 1 | id,family=poisson(log), verbose = FALSE) 

             beta <- as.vector(fit0$coeff$fixed[2:nitem])
             b <- as.vector(fit0$coeff$fixed[1]+fit0$coeff$random$id)
			 out$beta <- beta
			 out$b <- b
             out$mu <- fit0$coeff$fixed[1]
             out$sigma2 <- getVarCov(fit0)[1]
             return(out)
		   }
	 
		   if(is.null(offset))
	       {
			  roffset <- matrix(0,nrow=nsubject,ncol=p)
		   }
  	       else
		   {
			  roffset <- offset
		   }
	 
		   fit0 <- RaschMLE(ywork,p,nsubject,roffset)

         #########################################################################################
         # parameters depending on status
         #########################################################################################
       
     	   if(status==TRUE)
	       {
              beta <- fit0$beta
              b <- fit0$b

              muclus <- rep(0,maxn)
              sigmaclus <- rep(0,maxn)
              muclus[1] <- fit0$mu
              sigmaclus[1] <- fit0$sigma2

			  ncluster <- 1
              ss <- rep(1,nsubject)

	          if(homo.var==0)tauk2 <- 2.01

			  ccluster <- rep(0,maxn)
			  ccluster[1] <- nsubject
			  workv <- rep(0,maxn+1) 
			  wdp <- rep(0,maxn)
			  vdp <- rep(0,maxn)

			  foo <- .Fortran("dpweightsimbl",
                      n=as.integer(maxn),
                      ccluster=as.integer(ccluster),
                      alpha=as.double(alpha), 
                      workv=as.double(workv),
                      w=as.double(wdp),
                      v=as.double(vdp),
                      PACKAGE="DPpackage") 

			  wdp <- foo$w
			  vdp <- foo$v

   	       }
      	   if(status==FALSE)
	       {
			  alpha <- state$alpha
              beta <- state$beta
              b <- state$b
              ncluster <- state$ncluster
              ss <- state$ss
              muclus <- state$muclus   
              sigmaclus <- state$sigmaclus
              mub <- state$mub
              sigmab <- state$sigmab
              if(homo.var==0) tauk2 <- state$tauk2
              wdp <- state$wdp 
              vdp <- state$vdp
           }

         #########################################################################################
         # working space
         #########################################################################################

           seed1<-sample(1:29000,1)
           seed2<-sample(1:29000,1)
           seed<-c(seed1,seed2)

           betac <- rep(0,(p-1))
           iflagp <- rep(0,(p-1))
           xtx <- matrix(0,nrow=(p-1),ncol=(p-1))
           xty <- rep(0,(p-1))
           workmhp <- rep(0,(p-1)*p/2)
           workvp <- rep(0,p-1)

		   cstrt <- matrix(0,nrow=maxn,ncol=nsubject)
		   ccluster <- rep(0,maxn)
		   prob <- rep(0,maxn) 
		   workv <- rep(0,maxn+1)

           if(homo.var==0)
           {
              foo <- .Fortran("dpmraschp",
                            datastrm   =as.integer(datastrm),
  	 	                    imiss      =as.integer(imiss),
   	 	                    nmissing   =as.integer(nmissing),
							nsubject   =as.integer(nsubject),
	 	                    p          =as.integer(p),
  	 	                    y          =as.integer(ywork),
                            roffset    =as.double(roffset),
                            ngrid      =as.integer(ngrid),
                            grid       =as.double(grid),
  	 	                    maxn       =as.integer(maxn),  	 	
  	 	                    a0b0       =as.double(a0b0),
  	 	                    m0         =as.double(m0),
                            s0         =as.double(s0),
							prec       =as.double(prec),
  	 	                    sb         =as.double(sb),
							tau        =as.double(tau),
 		                    mcmc       =as.integer(mcmcvec),
							nsave      =as.integer(nsave),
                            acrate     =as.double(acrate),
 		                    cpo        =as.double(cpo),
							cpov       =as.double(cpov),
 		                    randsave   =as.double(randsave),
 		                    thetasave  =as.double(thetasave),
                            densave    =as.double(densave),
                            cdfsave    =as.double(cdfsave),
                            alpha      =as.double(alpha),		
 	 	                    b          =as.double(b),		
 		                    beta       =as.double(beta),	
                            ss         =as.integer(ss),
							ncluster   =as.integer(ncluster),
							mub        =as.double(mub),
 		                    sigmab     =as.double(sigmab),
                            tauk2      =as.double(tauk2),
                            wdp        =as.double(wdp),
                            vdp        =as.double(vdp),
						    muclus     =as.double(muclus),
							sigmaclus  =as.double(sigmaclus),
 		                    betac      =as.double(betac),		
 		                    workvp     =as.double(workvp),
                            workmhp    =as.double(workmhp), 
                            xtx        =as.double(xtx),
                            xty        =as.double(xty),
                            iflagp     =as.integer(iflagp),
                            cstrt      =as.integer(cstrt),
                            ccluster   =as.integer(ccluster),
                            prob       =as.double(prob),
                            workv      =as.double(workv),
                            seed       =as.integer(seed),
		                    PACKAGE    ="DPpackage")
           }
           else
           {
              foo <- .Fortran("dpmraschph",
                            datastrm   =as.integer(datastrm),
  	 	                    imiss      =as.integer(imiss),
   	 	                    nmissing   =as.integer(nmissing),
							nsubject   =as.integer(nsubject),
	 	                    p          =as.integer(p),
  	 	                    y          =as.integer(ywork),
                            roffset    =as.double(roffset),
                            ngrid      =as.integer(ngrid),
                            grid       =as.double(grid),
  	 	                    maxn       =as.integer(maxn),  	 	
  	 	                    a0b0       =as.double(a0b0),
  	 	                    m0         =as.double(m0),
                            s0         =as.double(s0),
							prec       =as.double(prec),
  	 	                    sb         =as.double(sb),
							tau        =as.double(tau),
 		                    mcmc       =as.integer(mcmcvec),
							nsave      =as.integer(nsave),
                            acrate     =as.double(acrate),
 		                    cpo        =as.double(cpo),
							cpov       =as.double(cpov),
 		                    randsave   =as.double(randsave),
 		                    thetasave  =as.double(thetasave),
                            densave    =as.double(densave),
                            cdfsave    =as.double(cdfsave),
                            alpha      =as.double(alpha),		
 	 	                    b          =as.double(b),		
 		                    beta       =as.double(beta),	
                            ss         =as.integer(ss),
							ncluster   =as.integer(ncluster),
							mub        =as.double(mub),
 		                    sigmab     =as.double(sigmab),
                            tauk2      =as.double(tauk2),
                            wdp        =as.double(wdp),
                            vdp        =as.double(vdp),
						    muclus     =as.double(muclus),
							sigmaclus  =as.double(sigmaclus),
 		                    betac      =as.double(betac),		
 		                    workvp     =as.double(workvp),
                            workmhp    =as.double(workmhp), 
                            xtx        =as.double(xtx),
                            xty        =as.double(xty),
                            iflagp     =as.integer(iflagp),
                            cstrt      =as.integer(cstrt),
                            ccluster   =as.integer(ccluster),
                            prob       =as.double(prob),
                            workv      =as.double(workv),
                            seed       =as.integer(seed),
		                    PACKAGE    ="DPpackage")
           }

         #########################################################################################
         # save state
         #########################################################################################

           hpdf <- function(x)
           {
                alpha <- 0.05
                vec <- x
                n <- length(x)         
                alow <- rep(0,2)
                aupp <- rep(0,2)
                a <-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                           alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                return(c(a$alow[1],a$aupp[1]))
           }
    
           pdf <- function(x)
           {
                alpha <- 0.05
                vec <- x
                n <- length(x)         
                alow<-rep(0,2)
                aupp<-rep(0,2)
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                          alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
              return(c(a$alow[2],a$aupp[2]))
           }

           model.name<-"Bayesian Rasch Poisson Count Model using a DPM of normals prior"		

           state <- list(alpha=foo$alpha,
	                     b=foo$b,
	                     beta=foo$beta,
                         ncluster=foo$ncluster,
                         ss=foo$ss,
						 muclus=foo$muclus,
                         sigmaclus=foo$sigmaclus,
	                     mub=foo$mub,
	                     sigmab=foo$sigmab,
                         tauk2=foo$tauk2,
                         wdp=foo$wdp,
                         vdp=foo$vdp)

           cpo <- matrix(foo$cpo,nrow=nsubject,ncol=p)
           randsave <- matrix(foo$randsave,nrow=nsave,ncol=nsubject+1)
           thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=p+6)
           densave <- matrix(foo$densave,nrow=nsave,ncol=ngrid)
           cdfsave <- matrix(foo$cdfsave,nrow=nsave,ncol=ngrid)

           dens.m <- apply(densave,2,mean)
           dens.l <- NULL
           dens.u <- NULL
           if(compute.band)
           {
              limm <- apply(densave, 2, hpdf)
              dens.l <- limm[1,]
              dens.u <- limm[2,]
           }

           cdf.m <- apply(cdfsave,2,mean)
           cdf.l <- NULL
           cdf.u <- NULL
           if(compute.band)
           {
              limm <- apply(cdfsave, 2, hpdf)
              cdf.l <- limm[1,]
              cdf.u <- limm[2,]
           }

           pnames <- paste("beta",2:p,sep="")
		   pnames <- c(pnames,"mu","sigma2","ncluster","alpha","mub","sigma2b","tauk2")
           colnames(thetasave) <- pnames
         
           qnames <- NULL
           for(i in 1:nsubject)
           {
               temp <- paste("theta(ID=",i,sep="")
               temp <- paste(temp,")",sep="")
               qnames <- c(qnames,temp)
           }
           qnames <- c(qnames,"prediction")

           dimnames(randsave) <- list(NULL,qnames)
         
           coeff <- apply(thetasave, 2, mean)
         
           save.state <- list(thetasave=thetasave,randsave=randsave,densave=densave,cdfsave=cdfsave)

         
		   acrate <- foo$acrate
         
	       z <- list(call=cl,
                     y=ywork,
                     modelname=model.name,
                     cpo=cpo,
                     prior=prior,
                     mcmc=mcmc, 
                     state=state,
                     save.state=save.state,
                     nrec=nrec,
                     nsubject=nsubject,
                     p=p,
                     acrate=acrate,
					 coefficients=coeff,
                     dens.m=dens.m,
                     dens.l=dens.l,
                     dens.u=dens.u,
                     cdf.m=cdf.m,
                     cdf.l=cdf.l,
                     cdf.u=cdf.u,
					 grid=grid,
                     alpharand=alpharand,
					 murand=murand,
                     sigmarand=sigmarand,
					 tauk2rand=tauk2rand,
                     cpov=foo$cpov)
                 
           cat("\n\n")
 	       class(z)<-c("DPMraschpoisson")
  	       return(z)
}



###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2009
### Last modification: 25-09-2009.
###



"print.DPMraschpoisson" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) 
    {
        cat("Posterior Inference of Parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)    
    }
    else cat("No coefficients\n")

    if(!is.null(x$acrate))
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   
   
    cat("\nNumber of subjects:",x$nsubject)
    cat("\nNumber of items:",x$p)
    cat("\n\n")
    invisible(x)
}

"summary.DPMraschpoisson"<-function(object, hpd=TRUE, ...) 
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

### Difficulty parameters

    dimen1 <- object$p-1

    if(dimen1>1)
    {
       mat <- thetasave[,1:dimen1]
    }
    else
    {
	   mat <- matrix(thetasave[,1:dimen1],ncol=1)
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

    names(coef.m) <- names(object$coefficients[1:dimen1])
    names(coef.sd) <- names(object$coefficients[1:dimen1])
    names(coef.se) <- names(object$coefficients[1:dimen1])
    names(coef.l) <- names(object$coefficients[1:dimen1])
    names(coef.u) <- names(object$coefficients[1:dimen1])

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
    
    ans <- c(object[c("call", "modelname")])

    ans$coefficients <- coef.table


### functionals

    mat<-thetasave[,(dimen1+1):(dimen1+2)]

    coef.p<-object$coefficients[(dimen1+1):(dimen1+2)]
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

    names(coef.m)<-names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.sd)<-names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.se)<-names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.l)<-names(object$coefficients[(dimen1+1):(dimen1+2)])
    names(coef.u)<-names(object$coefficients[(dimen1+1):(dimen1+2)])

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
    
    ans$functionals <-coef.table


### CPO
    ans$cpo<-object$cpo


### Precision parameter

    if(is.null(object$prior$a0))
    {
	   coef.p <- object$coefficients[(dimen1+3)] 
       mat <- matrix(thetasave[,(dimen1+3)],ncol=1)
    }
    else
    {
       coef.p <- object$coefficients[(dimen1+3:dimen1+4)] 
       mat <- thetasave[,(dimen1+3:dimen1+4)]
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

### Baseline Information


	dimen1 <- dimen1 + 4
    if((object$murand+object$sigmarand+object$tauk2rand)>0)
    {
       mat <- NULL
	   coef.p <- NULL
       if(object$murand==1)
   	   {
	       mat <- cbind(mat,thetasave[,dimen1+1]) 
	       coef.p <- c(coef.p,object$coefficients[dimen1+1])
       }
       if(object$sigmarand==1)
   	   {
	       mat <- cbind(mat,thetasave[,dimen1+2]) 
	       coef.p <- c(coef.p,object$coefficients[dimen1+2])
       }
       if(object$tauk2rand==1)
   	   {
	       mat <- cbind(mat,thetasave[,dimen1+3]) 
	       coef.p <- c(coef.p,object$coefficients[dimen1+3])
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

       ans$base<-coef.table
    }

    ans$nrec <- object$nrec
    ans$nsubject <- object$nsubject
    ans$p <- object$p
    ans$acrate <- object$acrate

    class(ans) <- "summaryDPMraschpoisson"
    return(ans)
}


"print.summaryDPMraschpoisson"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nDifficulty parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$functionals)) {
        cat("\nFunctionals of the Random Effects:\n")
        print.default(format(x$functionals, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

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

    if(!is.null(x$acrate))
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   

    cat("\nNumber of subjects:",x$nsubject)
    cat("\nNumber of items:",x$p,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.DPMraschpoisson"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPMraschpoisson"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients[1:(x$p+2)]
           n <- length(coef.p)
           pnames <- names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:n)
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
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

           if(x$alpharand==1)
           {
               title1 <- paste("Trace of","alpha",sep=" ")
               title2 <- paste("Density of","alpha",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+3]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+3],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }


           if(x$murand==1)
           {
               title1 <- paste("Trace of","mub",sep=" ")
               title2 <- paste("Density of","mub",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+4]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+4],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$sigmarand==1)
           {
               title1 <- paste("Trace of","sigma2b",sep=" ")
               title2 <- paste("Density of","sigma2b",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+5]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+5],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$tauk2rand==1)
           {
		      title1 <- paste("Trace of","tauk2",sep=" ")
              title2 <- paste("Density of","tauk2",sep=" ")       
		      plot(ts(x$save.state$thetasave[,x$p+6]),main=title1,xlab="MCMC scan",ylab=" ")
              fancydensplot1(x$save.state$thetasave[,x$p+6],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           title1 <- c("Density Estimate")
           title2 <- c("CDF Estimate")
           plot(x$grid,x$dens.m,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="theta")
           plot(x$grid,x$cdf.m,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")

        }
   
        else
        {
			coef.p <- x$coefficients
			n <- length(coef.p)
			pnames <- names(coef.p)
			poss <- 0 
			for(i in 1:n)
			{
				if(pnames[i]==param)poss <- i
			}
            if(poss==0 && param !="predictive")             
			{
				stop("This parameter is not present in the original model.\n")
			}
	    
			par(ask = ask)
			layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))

			if(param !="predictive")
			{
				par(ask = ask)
				layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
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
				title1<-c("Density Estimate")
				title2<-c("CDF Estimate")
				plot(x$grid,x$dens.m,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="theta")
				plot(x$grid,x$cdf.m,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")
            }                

        }
   }

}
