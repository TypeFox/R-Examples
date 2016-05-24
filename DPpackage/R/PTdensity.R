### PTdensity.R                   
### Fit a Mixture of Polya trees for density estimation
###
### Copyright: Alejandro Jara and Tim Hanson, 2006-2012.
###
### Last modification: 25-01-2010.
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
### The authors's contact information:
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
###      Tim Hanson
###      Department of Statistics
###      University of South Carolina
###      216 LeConte College, 
###      1523 Greene St
###      Columbia, SC 29208
###      Voice: 803-777-3853  URL  : http://www.stat.sc.edu/~hansont/
###      Fax  : 803-777-4048  Email: hansont@stat.sc.edu
###


PTdensity <- function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PTdensity")

PTdensity.default <- function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,data,na.action=na.fail)
{
         #########################################################################################
         # call parameters
         #########################################################################################
           cl <- match.call()
	       y <- na.action(as.matrix(y))	
	  
         #########################################################################################
         # data structure
         #########################################################################################
           nrec <- nrow(y)
           nvar <- ncol(y)
          
           if(nvar>1)
           {
              left <- rep(0,2)
	          right <- rep(0,2)
	    
			  left[1] <- min(y[,1])-0.5*sqrt(var(y[,1]))
	          right[1] <- max(y[,1])+0.5*sqrt(var(y[,1]))
	    
	          left[2] <- min(y[,2])-0.5*sqrt(var(y[,2]))
	          right[2] <- max(y[,2])+0.5*sqrt(var(y[,2]))
           } 
           else
           {
              left <- min(y)-0.5*sqrt(var(y))
              right <- max(y)+0.5*sqrt(var(y))
           }    

         #########################################################################################
         # prior information
         #########################################################################################

   		   jfr <- c(0,0)
           if(nvar==1)
           {
              tau <- c(-1,1)
              m0 <- 0
              S0 <- 1
           } 
           else
           {
			  nu0 <- -1
              tinv <- diag(1,nvar)
              m0 <- rep(0,nvar)
              S0 <- diag(1,nvar)
           }

           if(is.null(prior$mu))
           {
			  murand <- 1
              if(is.null(prior$m0))
              {
                 jfr[1] <- 1
              }
              else
              {
                 m0 <- prior$m0
                 S0 <- prior$S0
			  }
           }
           else
           {  
			  murand <- 0
              mu <- prior$mu
           }

           if(is.null(prior$sigma))
           {
			  sigmarand <- 1
              if(nvar==1)
              { 
                 if(is.null(prior$tau1))
                 {
                     jfr[2] <- 1
                 }
                 else
                 {
                     tau <- c(prior$tau1,prior$tau2)
                 }
              }
              else
              {
                 if(is.null(prior$nu0))
                 {
                    jfr[2] <- 1
                 }
                 else
                 {
                    nu0 <- prior$nu0
                    tinv <- prior$tinv
				 }   
              }
           }
           else
           {  
			  sigmarand <- 0
              sigma <- prior$sigma
		   }

  	       if(is.null(prior$a0))
  	       {
  	          ca <- -1
   	          cb <- -1 
  	          cpar <- prior$alpha
  	          crand <- 0
  	       }
           else
           {
              ca <- prior$a0
  	          cb <- prior$b0
  	          cpar <- 10
   	          crand <- 1
  	       }
  	       ab <- c(ca,cb)
  	 
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

           if(is.null(mcmc$tune3))         
           {
               tune3 <- 1.1
           }
           else
           {
               tune3 <- mcmc$tune3
           }


         #########################################################################################
         # output
         #########################################################################################
           acrate <- rep(0,3)

           if(nvar==1)
           {
               f <- rep(0,ngrid)
           }
           else
           {
               ngrid <- as.integer(sqrt(ngrid))
               f <- matrix(0,nrow=ngrid,ncol=ngrid)
               fun1 <- rep(0,ngrid)
               fun2 <- rep(0,ngrid)
           }
         
           thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+1)
           randsave <- matrix(0,nrow=nsave,ncol=nvar)
         
         #########################################################################################
         # parameters depending on status
         #########################################################################################
         
      	   if(status==TRUE)
	       {
	          if(nvar==1)
	          {
                 if(murand==1)
                 {
	                mu <- mean(y)
                 }
                 if(sigmarand==1)
                 {
	                sigma <- sqrt(var(y))
                 }
                 
	          }
	          else
	          {
                 if(murand==1)
                 {
					mu <- apply(y,2,mean)
                 }
                 if(sigmarand==1)
                 {
                    sigma <- matrix(var(y),nvar,nvar)
                 }
              }   
   	       }
	 
      	   if(status==FALSE)
	       {
	          cpar <- state$alpha
			  mu <- state$mu 
	          sigma <- state$sigma
	       }    

         #########################################################################################
         # working space
         #########################################################################################

           acrate <- rep(0,3)
	       cpo <- rep(0,nrec)
	       if(nvar==1)
		   {
              if(is.null(grid))
              {
	             grid <- seq(left,right,length=ngrid)
              }
              else
              {
                 grid <- as.matrix(grid)
                 ngrid <- nrow(grid)
                 grid <- as.vector(grid)
			  }
	       }        
	       else
	       {	
              if(is.null(grid))
              {
				 grid1<-seq(left[1],right[1],length=ngrid)  
                 grid2<-seq(left[2],right[2],length=ngrid)  
              }
              else
              {
                 grid <- as.matrix(grid)
                 ngrid <- nrow(grid)
                 grid1 <- grid[,1]
                 grid2 <- grid[,2]
			  }
           }

	       seed <- c(sample(1:29000,1),sample(1:29000,1))

         #########################################################################################
         # calling the fortran code
         #########################################################################################

           if(nvar==1)
           {
              if(is.null(prior$M))
              {
				  whicho <- rep(0,nrec)
    	          whichn <- rep(0,nrec)
				  foo <- .Fortran("ptdensityu",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						y          =as.double(y),
						ab         =as.double(ab),
						murand     =as.integer(murand),
						sigmarand  =as.integer(sigmarand),
						jfr        =as.integer(jfr),
						m0         =as.double(m0),
						s0         =as.double(S0),  
						tau        =as.double(tau),
						mcmcvec    =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						tune1      =as.double(tune1),
						tune2      =as.double(tune2),
						tune3      =as.double(tune3),
						acrate     =as.double(acrate),
						f          =as.double(f),
						thetasave  =as.double(thetasave),		
						cpo        =as.double(cpo),		
						cpar       =as.double(cpar),		
						mu         =as.double(mu),		
						sigma      =as.double(sigma),		
						grid       =as.double(grid),		
						seed       =as.integer(seed),
						whicho     =as.integer(whicho),
						whichn     =as.integer(whichn),
						PACKAGE    ="DPpackage")
			}	
            else
            {
				  nlevel <- prior$M
                  ninter <- 2**nlevel
                  assign <- matrix(0,nrow=nrec,ncol=nlevel)
				  accums <- matrix(0,nrow=nlevel,ncol=ninter)
                  counter <- matrix(0,nrow=nlevel,ncol=ninter)
                  endp <- rep(0,ninter-1)
				  intpn <- rep(0,nrec)
				  intpo <- rep(0,nrec)
				  prob <- rep(0,ninter)
                  rvecs <- matrix(0,nrow=nlevel,ncol=ninter)
 
				  foo <- .Fortran("ptdensityup",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						y          =as.double(y),
						ab         =as.double(ab),
						murand     =as.integer(murand),
						sigmarand  =as.integer(sigmarand),
						jfr        =as.integer(jfr),
						m0         =as.double(m0),
						s0         =as.double(S0),  
						tau        =as.double(tau),
						nlevel     =as.integer(nlevel),
						ninter     =as.integer(ninter),
						mcmcvec    =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						tune1      =as.double(tune1),
						tune2      =as.double(tune2),
						tune3      =as.double(tune3),
						acrate     =as.double(acrate),
						f          =as.double(f),
						thetasave  =as.double(thetasave),		
						cpo        =as.double(cpo),		
						cpar       =as.double(cpar),		
						mu         =as.double(mu),		
						sigma      =as.double(sigma),		
						grid       =as.double(grid),		
						intpn     =as.integer(intpn),		
						intpo     =as.integer(intpo),		
						accums    =as.double(accums),
						assign    =as.integer(assign),
						counter   =as.integer(counter),
						endp      =as.double(endp),
						prob      =as.double(prob),
						rvecs     =as.double(rvecs),
						seed      =as.integer(seed),
						PACKAGE    ="DPpackage")

			  }
	       }   
           else
           {
				iflag <- rep(0,nvar)
				limw <- rep(0,nvar)
				linf <- rep(0,nvar)
				lsup <- rep(0,nvar)
				muc <- rep(0,nvar)
				narea <- 2**nvar
				mass <- rep(0,narea)
				massi <- rep(0,narea)
				parti <- rep(0,nvar)
				pattern <- rep(0,nvar)
				patterns <- rep(0,nvar)
				propv <- matrix(0,nrow=nvar,ncol=nvar) 
				propv1 <- matrix(0,nrow=nvar,ncol=nvar) 
				propv2 <- matrix(0,nrow=nvar,ncol=nvar) 
				s <- matrix(0,nrow=nvar,ncol=nvar)         
				sigmac <- matrix(0,nrow=nvar,ncol=nvar)         
				sigmainv <- matrix(0,nrow=nvar,ncol=nvar)         
				sigmainvc <- matrix(0,nrow=nvar,ncol=nvar)         
				vv <- rep(0,nvar)         
				whicho <- rep(0,nrec)
				whichn <- rep(0,nrec)
				workh1 <- rep(0,nvar*(nvar+1)/2)
				workh2 <- rep(0,nvar*(nvar+1)/2)
				workmh <- rep(0,nvar*(nvar+1)/2)
				workm1 <- matrix(0,nrow=nvar,ncol=nvar)         
				workm2 <- matrix(0,nrow=nvar,ncol=nvar)         
				ybar <- rep(0,nvar) 
				z <- matrix(0,nrow=nrec,ncol=nvar)
				zc <- matrix(0,nrow=nrec,ncol=nvar)
				zwork <- rep(0,nvar)         

				s0 <- solve(S0)
         
				if(is.null(prior$M))
				{
					foo <- .Fortran("ptmdensity",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						ab         =as.double(ab),
						murand     =as.integer(murand),
						sigmarand  =as.integer(sigmarand),
						jfr        =as.integer(jfr),
						m0         =as.double(m0),
						s0         =as.double(s0),  
						nu0        =as.integer(nu0),
						tinv	   =as.double(tinv),
						mcmcvec    =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						tune1      =as.double(tune1),
						tune2      =as.double(tune2),
						tune3      =as.double(tune3),
						acrate     =as.double(acrate),
						cpo        =as.double(cpo),		
						f          =as.double(f),
						randsave   =as.double(randsave),		
						thetasave  =as.double(thetasave),		
						cpar       =as.double(cpar),		
						mu         =as.double(mu),		
						sigma      =as.double(sigma),		
						grid1      =as.double(grid1),
						grid2      =as.double(grid2),
						iflag      =as.integer(iflag),
						whicho     =as.integer(whicho),
						whichn     =as.integer(whichn),
						limw       =as.double(limw),
						linf       =as.double(linf),
						lsup       =as.double(lsup),
						narea      =as.integer(narea),
						mass       =as.double(mass),
						massi      =as.integer(massi),
						parti      =as.integer(parti),
						pattern    =as.integer(pattern),
						patterns   =as.integer(patterns),
						s          =as.double(s),
						sigmainv   =as.double(sigmainv),
						sigmainvc  =as.double(sigmainvc),
						ybar       =as.double(ybar),
						z          =as.double(z),
						zc         =as.double(zc),
						zwork      =as.double(zwork),
						vv         =as.double(vv),
						workmh     =as.double(workmh),
						workh1     =as.double(workh1),
						workh2     =as.double(workh2),
						workm1     =as.double(workm1),
						workm2     =as.double(workm2),
						muc        =as.double(muc),
						sigmac     =as.double(sigmac),
						propv      =as.double(propv),
						propv1     =as.double(propv1),
						propv2     =as.double(propv2),
						seed       =as.integer(seed),
						PACKAGE    ="DPpackage")
				}
				else
				{
					foo <- .Fortran("ptmdensityp",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						ab         =as.double(ab),
                        nlevel     =as.integer(prior$M),
						murand     =as.integer(murand),
						sigmarand  =as.integer(sigmarand),
						jfr        =as.integer(jfr),
						m0         =as.double(m0),
						s0         =as.double(s0),  
						nu0        =as.integer(nu0),
						tinv	   =as.double(tinv),
						mcmcvec    =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						tune1      =as.double(tune1),
						tune2      =as.double(tune2),
						tune3      =as.double(tune3),
						acrate     =as.double(acrate),
						cpo        =as.double(cpo),		
						f          =as.double(f),
						randsave   =as.double(randsave),		
						thetasave  =as.double(thetasave),		
						cpar       =as.double(cpar),		
						mu         =as.double(mu),		
						sigma      =as.double(sigma),		
						grid1      =as.double(grid1),
						grid2      =as.double(grid2),
						iflag      =as.integer(iflag),
						whicho     =as.integer(whicho),
						whichn     =as.integer(whichn),
						limw       =as.double(limw),
						linf       =as.double(linf),
						lsup       =as.double(lsup),
						narea      =as.integer(narea),
						mass       =as.double(mass),
						massi      =as.integer(massi),
						parti      =as.integer(parti),
						pattern    =as.integer(pattern),
						patterns   =as.integer(patterns),
						s          =as.double(s),
						sigmainv   =as.double(sigmainv),
						sigmainvc  =as.double(sigmainvc),
						ybar       =as.double(ybar),
						z          =as.double(z),
						zc         =as.double(zc),
						zwork      =as.double(zwork),
						vv         =as.double(vv),
						workmh     =as.double(workmh),
						workh1     =as.double(workh1),
						workh2     =as.double(workh2),
						workm1     =as.double(workm1),
						workm2     =as.double(workm2),
						muc        =as.double(muc),
						sigmac     =as.double(sigmac),
						propv      =as.double(propv),
						propv1     =as.double(propv1),
						propv2     =as.double(propv2),
						seed       =as.integer(seed),
						PACKAGE    ="DPpackage")
				}
           }

         #########################################################################################
         # save state
         #########################################################################################
           model.name <- "Bayesian Density Estimation Using MPT"		
                
           varnames<-colnames(y)
           if(is.null(varnames))
           {
               varnames<-all.vars(cl)[1:nvar]
		   }
		
           state <- list(alpha=foo$cpar,
	                     mu=foo$mu,
	                     sigma=matrix(foo$sigma,nrow=nvar,ncol=nvar)
                         )
         
		   thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=(nvar+nvar*(nvar+1)/2+1))
           if(nvar>1)
           {
              randsave <- matrix(foo$randsave,nrow=nsave,ncol=nvar)
              colnames(randsave) <- varnames
           }   

           coeff<-apply(thetasave,2,mean) 
         
           pnames1<-NULL
	       for(i in 1:nvar)
	       {
	           pnames1<-c(pnames1,paste("mu",varnames[i],sep=":"))
		   }
           pnames2<-NULL
      	   for(i in 1:nvar)
	       {
	           for(j in i:nvar)
	           {
	               if(i==j)
	               {
	                  tmp <- varnames[i]
				   }
				   else
	               {
	                  tmp <- paste(varnames[i],varnames[j],sep="-")
	               }   
				   pnames2 <- c(pnames2,paste("sigma",tmp,sep=":"))
			   }	
		   }
	 
           names(coeff) <- c(pnames1,pnames2,"alpha")
           colnames(thetasave) <- c(pnames1,pnames2,"alpha")
           save.state <- list(thetasave=thetasave,randsave=randsave)
         
		   if(crand==0)
           {
              acrate<-foo$acrate[1:2]
           }
           else
           {
            acrate<-foo$acrate
           }

           x1<-NULL
           x2<-NULL
           dens<-NULL
         
         if(nvar==1)
         {
            x1<-foo$grid
            dens<-foo$f
            f<-foo$f
            grid1<-foo$grid
            grid2<-NULL
            fun1<-foo$f
            fun2<-NULL
         }
         else
         {
            x1 <- grid1
            x2 <- grid2
            dens <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
            f <- matrix(foo$f,nrow=ngrid,ncol=ngrid)

            dist1 <- grid2[2]-grid2[1] 
			dist2 <- grid1[2]-grid1[1] 
            fun1 <- (dist1/2)*(dens[,1]+dens[,ngrid]+2*apply(dens[,2:(ngrid-1)],1,sum))
            fun2 <- (dist2/2)*(dens[1,]+dens[ngrid,]+2*apply(dens[2:(ngrid-1),],2,sum))
         }   

	 z<-list(call=cl,y=y,varnames=varnames,modelname=model.name,cpo=foo$cpo,
                 prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
                 nvar=nvar,crand=crand,coefficients=coeff,f=f,grid1=grid1,grid2=grid2,
                 fun1=fun1,fun2=fun2,x1=x1,x2=x2,dens=dens,acrate=acrate)
                 
         cat("\n\n")
 	 class(z)<-"PTdensity"
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 28-11-2006.
###


"print.PTdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Variables:",x$nvar,"\n")        
    cat("\n\n")
    invisible(x)
}


"summary.PTdensity"<-function(object, hpd=TRUE, ...) 
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
    nvar<-object$nvar

    if(object$crand==0)
    {
       dimen<-(nvar+nvar*(nvar+1)/2)
       mat<-matrix(thetasave[,1:dimen],ncol=2) 
       
    }
    else
    {
       dimen<-(nvar+nvar*(nvar+1)/2+1)
       mat<-thetasave[,1:dimen]
    }

    coef.p<-object$coefficients[1:dimen]
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

    names(coef.m)<-names(object$coefficients[1:dimen])
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.se)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])

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

    ans$coefficients<-coef.table


### CPO
    ans$cpo<-object$cpo

    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec
    ans$nvar<-object$nvar

    class(ans) <- "summaryPTdensity"
    return(ans)
}


"print.summaryPTdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nBaseline parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Variables:",x$nvar,"\n")            
    cat("\n\n")
    invisible(x)
}


"plot.PTdensity"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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


   if(is(x, "PTdensity"))
   {

      if(output=="density")
      {

      # Density estimation
	
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))

        if(x$nvar==1)
        {
           title1<-paste("Density of",x$varnames[1],sep=' ')
	       aa<-hist(x$y[,1],plot=F,)
		   maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
		   miny<-min(x$y[,1])
	       maxy<-max(x$y[,1])
	       deltay<-(maxy-miny)*0.2
	       miny<-miny-deltay
	       maxy<-maxy+deltay
	      
	        hist(x$y[,1],probability=T,xlim=c(min(x$grid1),max(x$grid1)),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")
           lines(x$x1,x$dens,lwd=2)
        }

        if(x$nvar==2)
        {
           title1<-paste("Density of",x$varnames[1],sep=' ')
	       aa<-hist(x$y[,1],plot=F,)
 	       maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
		   miny<-min(x$y[,1])
	       maxy<-max(x$y[,1])
	       deltay<-(maxy-miny)*0.2
	       miny<-miny-deltay
	       maxy<-maxy+deltay
	      
	       hist(x$y[,1],probability=T,xlim=c(min(x$grid1),max(x$grid1)),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")
           lines(x$grid1,x$fun1,lwd=2)

           title1<-paste("Density of",x$varnames[2],sep=' ')
		   aa<-hist(x$y[,2],plot=F,)
 	       maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
		   miny<-min(x$y[,2])
	       maxy<-max(x$y[,2])
	       deltay<-(maxy-miny)*0.2
	       miny<-miny-deltay
	       maxy<-maxy+deltay
	      
		   hist(x$y[,2],probability=T,xlim=c(min(x$grid2),max(x$grid2)),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")
           lines(x$grid2,x$fun2,lwd=2)

           varsn<-paste(x$varnames[1],x$varnames[2],sep="-")
	       title1<-paste("Predictive Density of ",varsn,sep='')

           xx<-matrix(x$grid1,ncol=1)
           yy<-matrix(x$grid2,ncol=1)
           z<-x$f 
           colnames(xx)<-x$varnames[1]
           colnames(yy)<-x$varnames[2]

	       contour(xx,yy,z,main=title1,xlab=x$varnames[1],ylab=x$varnames[2])
	       persp(xx,yy,z,xlab=x$varnames[1],ylab=x$varnames[2],zlab="density",theta=-30,phi=15,expand = 0.9, ltheta = 120,main=title1)
        }

        if(x$nvar>2)
        {
           for(i in 1:x$nvar)
    	   {
	           title1<-paste("Density of",x$varnames[i],sep=' ')
	    
	           aa<-hist(x$y[,i],plot=F)
	           maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
	           miny<-min(x$y[,i])-0.5*sqrt(var(x$y[,i]))
	           maxy<-max(x$y[,i])+0.5*sqrt(var(x$y[,i]))
	           deltay<-(maxy-miny)*0.2
	           miny<-miny-deltay
	           maxy<-maxy+deltay
	           hist(x$y[,i],probability=T,xlim=c(miny,maxy),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")
               lines(density(x$save.state$randsave[,i]),lwd=2)
	   }
	
           for(i in 1:(x$nvar-1))
		   {

               vectmp<-x$y[,i]
               xlim<-c(min(vectmp)-0.5*sqrt(var(vectmp)),max(vectmp)+0.125*sqrt(var(vectmp)))
	   
			   for(j in (i+1):x$nvar)
			   {
                   vectmp<-x$y[,j]
                   ylim<-c(min(vectmp)-0.5*sqrt(var(vectmp)),max(vectmp)+0.125*sqrt(var(vectmp)))

				   varsn<-paste(x$varnames[i],x$varnames[j],sep="-")
				   title1<-paste("Predictive Density of ",varsn,sep='')
                   xx<-x$save.state$randsave[,i]	    
                   yy<-x$save.state$randsave[,j]	    
	               est<-bivk(xx,yy,n=200)
	               contour(est,xlim=xlim,ylim=ylim,main=title1,xlab=x$varnames[i],ylab=x$varnames[j])
	               persp(est,theta=-30,phi=15,expand = 0.9, ltheta = 120,main=title1,
	                 xlab=x$varnames[i],ylab=x$varnames[j],zlab="density")
	   	}
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
               fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
            fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
        }


      }	
   }
}





