### LDDPraschpoisson.R                   
### Fit a Rasch Poisson model with a Linear Dependent DP prior
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


LDDPraschpoisson <- function(formula,prior,mcmc,offset=NULL,state,status,
			                 grid=seq(-10,10,length=1000),zpred,data=sys.frame(sys.parent()),compute.band=FALSE)
UseMethod("LDDPraschpoisson")

LDDPraschpoisson.default <-
function(formula,
         prior,
         mcmc,
         offset=NULL,         
         state,
         status,
         grid=seq(-10,10,length=1000),
         zpred,
         data=sys.frame(sys.parent()),
         compute.band=FALSE)
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

		   y <- model.response(mf,"numeric")

         #########################################################################################
         # data structure
         #########################################################################################
		   nsubject <- nrow(y)
	       p <- ncol(y)
           ywork <- y

		   z <- as.matrix(model.matrix(formula))
		   q <- ncol(z)

           if(ncol(zpred) != q)
           { 
			  stop("The design matrix for prediction must have the same number of columns than the data desing matrix.\n")     
           }

         #########################################################################################
         # prediction
         #########################################################################################
           npred <- nrow(zpred)
           ngrid <- length(grid)

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

           tau1 <- prior$tau1
           if(tau1 < 0)
           { 
			  stop("The parameter of the Gamma prior for the kernel variance must be possitive.\n")     
           }

           if(is.null(prior$tau2))
           {
              taus1 <- prior$taus1
              taus2 <- prior$taus2
              if(taus1 < 0 || taus2 < 0)
              { 
                  stop("The parameters of the Gamma prior for the gamma centering distribution must be possitive.\n")     
              }
			  tau2 <- 2.01
              tau2rand <- 1
           } 
           else
           {
			  taus1 <- -1
              taus2 <- -1
              tau2 <- prior$tau2
              tau2rand <- 0
           }

            b0 <- prior$beta0
            prec1 <- solve(prior$Sbeta0)
            sb <- prec1%*%b0

            if(length(b0)!=(p-1))
            { 
			   stop("Error in the dimension of the mean of the normal prior for the difficulty parameters.\n")     
            }

            if(nrow(prec1)!=(p-1) || ncol(prec1)!=(p-1))
			{ 
			   stop("Error in the dimension of the covariance of the normal prior for the difficulty parameters.\n")     
            }

			if(is.null(prior$mu0))
			{
				mu0 <- rep(0,q)
                mu <- prior$mub
				prec2 <- diag(1,q)
				murand <- 0
			}
			else
			{
				mu0 <- prior$mu0
				prec2 <- solve(prior$S0)
				mu <- rep(0,q)
				murand <- 1
			}
            smu <- prec2%*%mu0

            if(length(mu0)!=q)
            { 
			   stop("Error in the dimension of the prior mean of the mean of the normal centering distribution.\n")     
            }

            if(length(mu)!=q)
            { 
			   stop("Error in the dimension of the mean of the normal centering distribution.\n")     
            }

            if(nrow(prec2)!=q || ncol(prec2)!=q)
			{ 
			   stop("Error in the dimension of the variance of the mean of the normal centering distribution.\n")     
            }


            if(is.null(prior$nu))
            {
				nu <- -1
				tinv <- diag(1,q)
			    sigma <- prior$sb
				sigmarand <- 0
			}
			else
			{
				nu <- prior$nu
				tinv <- prior$psiinv
				sigma <- diag(1,q)
				sigmarand<-1
			}

            if(nrow(tinv)!=q || ncol(tinv)!=q)
			{ 
			   stop("Error in the dimension of the Wishart prior for the variance of the normal centering distribution.\n")     
            }

            if(nrow(sigma)!=q || ncol(sigma)!=q)
			{ 
			   stop("Error in the dimension of the normal cenetering covariance matrix.\n")     
            }

            a0b0 <- c(a0b0,tau1,taus1,taus2,nu)

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
		   acrate <- rep(0,2)
		   cpo <- matrix(0,nrow=nsubject,ncol=p)
		   denspm <- matrix(0,nrow=npred,ncol=ngrid)		
		   randsave <- matrix(0,nrow=nsave,ncol=nsubject+npred)

           thetasave <- matrix(0,nrow=nsave,ncol=p+q+(q*(q+1)/2)+2)
           densave <- matrix(0,nrow=nsave,ncol=npred*ngrid)

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

			  ncluster <- 1
			  ss <- rep(1,nsubject)

              alphaclus <- matrix(0,nrow=nsubject+100,ncol=q)
              sigmaclus <- rep(0,nsubject+100)
              alphaclus[1,1] <- fit0$mu
              sigmaclus[1] <- fit0$sigma2

   	       }
      	   if(status==FALSE)
	       {
              beta <- state$beta
              b <- state$b

			  alpha <- state$alpha
              ncluster <- state$ncluster
              ss <- state$ss
              alphaclus <- state$alphaclus   
              sigmaclus <- state$sigmaclus

              mu <- state$mub
              sigma <- state$sb
              tau2 <- state$tau2
           }

         #########################################################################################
         # working space
         #########################################################################################

           seed1 <- sample(1:29000,1)
           seed2 <- sample(1:29000,1)
           seed <- c(seed1,seed2)

		   iflagp <- rep(0,(p-1))
           betac <- rep(0,(p-1))
           xtx <- matrix(0,nrow=(p-1),ncol=(p-1))
           xty <- rep(0,(p-1))
           workmhp1	<- rep(0,(p-1)*p/2)
           workvp1 <- rep(0,p-1)

		   cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)
		   ccluster <- rep(0,nsubject)
		   iflagq <- rep(0,q)
		   alphawork <- rep(0,q)
           densw <- matrix(0,nrow=npred,ncol=ngrid)
		   prob <- rep(0,nsubject+100) 
		   quadf <- matrix(0,nrow=q,ncol=q)
		   sigmainv <- matrix(0,nrow=q,ncol=q)
		   workmhq1 <- rep(0,q*(q+1)/2)
		   workmhq2 <- rep(0,q*(q+1)/2)
		   workvq1 <- rep(0,q)
		   workvq2 <- rep(0,q)
		   ztz <- matrix(0,nrow=q,ncol=q)
		   zty <- rep(0,q)

         #########################################################################################
         # calling the Fortran code
         #########################################################################################

           foo <- .Fortran("lddpraschpoi",
							ngrid		=as.integer(ngrid),
							npred		=as.integer(npred),
							nsubject	=as.integer(nsubject),
							p			=as.integer(p),
							q			=as.integer(q),
							y			=as.integer(y),
							roffset		=as.double(roffset),
							grid		=as.double(grid),
							z			=as.double(z),
							zpred		=as.double(zpred),
							murand		=as.integer(murand),
							a0b0		=as.double(a0b0),
							b0			=as.double(b0),
							prec1		=as.double(prec1),
							sb			=as.double(sb),
							mu0			=as.double(mu0),
							prec2		=as.double(prec2),	
							smu			=as.double(smu),
							tinv		=as.double(tinv),
							acrate		=as.double(acrate),
							cpo			=as.double(cpo),
							denspm		=as.double(denspm),
							randsave	=as.double(randsave),
							thetasave	=as.double(thetasave),
							densave		=as.double(densave),
							ncluster	=as.integer(ncluster),
							ss			=as.integer(ss),
							beta		=as.double(beta),
							b			=as.double(b),
							alphaclus	=as.double(alphaclus),
							sigmaclus	=as.double(sigmaclus),
							alpha    	=as.double(alpha),
							mu			=as.double(mu),
							sigma		=as.double(sigma),
							tau2		=as.double(tau2),
							mcmc		=as.integer(mcmcvec),	
							nsave		=as.integer(nsave),
							iflagp		=as.integer(iflagp),
							betac		=as.double(betac),
							xtx			=as.double(xtx),
							xty			=as.double(xty),
							workmhp1	=as.double(workmhp1),
							workvp1		=as.double(workvp1),
							cstrt		=as.integer(cstrt),
							ccluster	=as.integer(ccluster),
							iflagq		=as.integer(iflagq),
							alphawork	=as.double(alphawork),
							densw		=as.double(densw),
							prob		=as.double(prob),
							quadf		=as.double(quadf),
							sigmainv	=as.double(sigmainv),
							workmhq1	=as.double(workmhq1),
							workmhq2	=as.double(workmhq2),
							workvq1		=as.double(workvq1),
							workvq2		=as.double(workvq2),
							ztz			=as.double(ztz),
							zty			=as.double(zty),
                            seed        =as.integer(seed),
		                    PACKAGE     ="DPpackage")


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

           model.name<-"Bayesian Semiparametric Rasch Poisson Model using a LDDP prior"		

           state <- list(alpha=foo$alpha,
	                     b=foo$b,
	                     beta=foo$beta,
                         ncluster=foo$ncluster,
                         ss=foo$ss,
						 alphaclus=matrix(foo$alphaclus,nrow=nsubject+100,ncol=q),
                         sigmaclus=foo$sigmaclus,
	                     mub=foo$mu,
	                     sb=matrix(foo$sigma,nrow=q,ncol=q),
                         tau2=foo$tau2)

           cpo <- matrix(foo$cpo,nrow=nsubject,ncol=p)
           randsave <- matrix(foo$randsave,nrow=nsave,ncol=nsubject+npred)
           thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=p+q+(q*(q+1)/2)+2)
           densave <- matrix(foo$densave,nrow=nsave,ncol=npred*ngrid)

           dens.m <- matrix(foo$denspm,nrow=npred,ncol=ngrid)		
           dens.l <- NULL
           dens.u <- NULL
           if(compute.band)
           {
              limm <- apply(densave, 2, hpdf)
              dens.l <- limm[1,]
              dens.u <- limm[2,]
			  dens.l <- matrix(dens.l,nrow=npred,ncol=ngrid,byrow=TRUE)
			  dens.u <- matrix(dens.u,nrow=npred,ncol=ngrid,byrow=TRUE)
           }

           pnames <- paste("beta",2:p,sep="")
		   pnames <- c(pnames,"tau2",paste("mub",seq(1,q),sep="-"))
           for(i in 1:q)
           {
			   for(j in i:q)
			   {
                   pnames <- c(pnames,paste("sb",i,j,sep=""))
			   }
		   }
		   pnames <- c(pnames,"ncluster","alpha")

           colnames(thetasave) <- pnames
         
           qnames <- NULL
           for(i in 1:nsubject)
           {
               temp <- paste("theta(ID=",i,sep="")
               temp <- paste(temp,")",sep="")
               qnames <- c(qnames,temp)
           }
           qnames <- c(qnames,paste("prediction",seq(1,npred),sep="-"))

           dimnames(randsave) <- list(NULL,qnames)
         
           coeff <- apply(thetasave, 2, mean)
         
           save.state <- list(thetasave=thetasave,
							  randsave=randsave,
							  densave=densave)
         
		   acrate <- foo$acrate
         
	       z <- list(call=cl,
                     y=ywork,
                     modelname=model.name,
                     cpo=cpo,
                     prior=prior,
                     mcmc=mcmc, 
                     state=state,
                     save.state=save.state,
                     nsubject=nsubject,
                     p=p,
					 q=q,
					 npred=npred,
					 zpred=zpred,
                     acrate=acrate,
					 coefficients=coeff,
                     dens.m=dens.m,
                     dens.l=dens.l,
                     dens.u=dens.u,
					 grid=grid,
                     alpharand=alpharand,
					 murand=murand,
                     sigmarand=sigmarand,
					 tau2rand=tau2rand,
					 compute.band=compute.band)
                 
           cat("\n\n")
 	       class(z) <- c("LDDPrasch")
  	       return(z)
}


