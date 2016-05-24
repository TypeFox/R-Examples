### DPdensity.R                   
### Fit a linear Dirichlet Process mixture of normal model for
### density estimation
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 05-10-2009.
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

DPdensity<-function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,method="neal",data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPdensity")

DPdensity.default <- function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,method="neal",data,na.action=na.fail)
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
           
           left <- rep(0,2)
           right <- rep(0,2)
         
           if(nvar==1)
           {
				left[1] <- min(y)-0.5*sqrt(var(y))
				right[1] <- max(y)+0.5*sqrt(var(y))
           }
           else
           {
				left[1] <- min(y[,1])-0.5*sqrt(var(y[,1]))
				right[1] <- max(y[,1])+0.5*sqrt(var(y[,1]))

				left[2] <- min(y[,2])-0.5*sqrt(var(y[,2]))
				right[2] <- max(y[,2])+0.5*sqrt(var(y[,2]))
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
         # output
         #########################################################################################

           if(is.null(grid))
           {
              if(nvar>1) ngrid <- as.integer(sqrt(ngrid))
			  grid1 <- seq(left[1],right[1],length=ngrid)  
			  grid2 <- seq(left[2],right[2],length=ngrid)  

		   }
           else
           {
			  grid <- as.matrix(grid)
			  ngrid <- nrow(grid)
              if(nvar==1)
              {
			     grid1 <- grid[,1]
				 grid2 <- rep(0,ngrid)
              }
              else
              {
				grid1 <- grid[,1]
				grid2 <- grid[,2]
			  }
		   }
		   f <- matrix(0,nrow=ngrid,ncol=ngrid)
           fun1 <- rep(0,ngrid)
           fun2 <- rep(0,ngrid)
           thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
           randsave <- matrix(0,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
         

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
         # working space
         #########################################################################################
           ccluster <- rep(0,nrec)
           cpo <- rep(0,nrec) 
           iflag <- rep(0,nvar) 
           muwork <- rep(0,nvar) 
           muwork2 <- rep(0,nvar) 
           prob <- rep(0,(nrec+100))
           s1 <- matrix(0,nvar,nvar)
           seed1 <- sample(1:29000,1)
           seed2 <- sample(1:29000,1)
           seed3 <- sample(1:29000,1)
           seed <- c(seed1,seed2,seed2)
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
           workcpo <- rep(0,nrec)

         #########################################################################################
         # calling the fortran code
         #########################################################################################
         
           if(method=="no-gaps")
           {
             if(nvar<=2)
             {
                foo <- .Fortran("bivspdeng",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						a0b0       =as.double(a0b0),
						k0         =as.double(k0),
						nuvec      =as.integer(nuvec),
						m1rand     =as.integer(m1rand),
						s2inv      =as.double(s2inv),
						s2invm2    =as.double(s2invm2),
						psiinv2    =as.double(psiinv2),
						tau        =as.double(tau),
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						cpo        =as.double(cpo),
						f          =as.double(f),
						fun1       =as.double(fun1), 		
						fun2       =as.double(fun2), 		 		
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
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
						grid1      =as.double(grid1),
						grid2      =as.double(grid2), 		
						iflag      =as.integer(iflag),
						muwork     =as.double(muwork),
						muwork2    =as.double(muwork2), 		
						prob       =as.double(prob),
						seed       =as.integer(seed),
						sigmawork  =as.double(sigmawork),
						sigmawork2 =as.double(sigmawork2), 		
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
						workcpo    =as.double(workcpo),
						PACKAGE    ="DPpackage")
			  }
			  else
			  {         
                foo <- .Fortran("spdeng",
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						a0b0       =as.double(a0b0),
						k0         =as.double(k0),
						nuvec      =as.integer(nuvec),
						m1rand     =as.integer(m1rand),
						s2inv      =as.double(s2inv),
						s2invm2    =as.double(s2invm2),
						psiinv2    =as.double(psiinv2),
						tau        =as.double(tau),
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
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
						iflag      =as.integer(iflag),
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
                        workcpo    =as.double(workcpo),
						PACKAGE    ="DPpackage")
               }		
           }		

           if(method=="neal")
           {
               if(nvar<=2)
               {
                foo <- .Fortran("bivspdenn",
						ngrid      =as.integer(ngrid),
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						a0b0       =as.double(a0b0),
						k0         =as.double(k0),
						nuvec      =as.integer(nuvec),
						m1rand     =as.integer(m1rand),
						s2inv      =as.double(s2inv),
						s2invm2    =as.double(s2invm2),
						psiinv2    =as.double(psiinv2),
						tau        =as.double(tau),
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						cpo        =as.double(cpo),
						f          =as.double(f),
						fun1       =as.double(fun1), 		
						fun2       =as.double(fun2), 		 		
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
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
						grid1      =as.double(grid1),
						grid2      =as.double(grid2), 		
						iflag      =as.integer(iflag),
						muwork     =as.double(muwork),
						muwork2    =as.double(muwork2), 		
						prob       =as.double(prob),
						seed       =as.integer(seed),
						sigmawork  =as.double(sigmawork),
						sigmawork2 =as.double(sigmawork2), 		
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
						workcpo    =as.double(workcpo),
						PACKAGE    ="DPpackage")
             }
             else
             {
                foo <- .Fortran("spdenn",
						nrec       =as.integer(nrec),
						nvar       =as.integer(nvar),
						y          =as.double(y),
						a0b0       =as.double(a0b0),
						k0         =as.double(k0),
						nuvec      =as.integer(nuvec),
						m1rand     =as.integer(m1rand),
						s2inv      =as.double(s2inv),
						s2invm2    =as.double(s2invm2),
						psiinv2    =as.double(psiinv2),
						tau        =as.double(tau),
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
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
						iflag      =as.integer(iflag),
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
						workcpo    =as.double(workcpo),
						PACKAGE    ="DPpackage")
               }
           }		

         #########################################################################################
         # save state
         #########################################################################################

           model.name <- "DPM of normals model for density estimation"		
         
		   varnames <- colnames(y)
	
           state <- list(
						alpha=foo$alpha,
						m1=matrix(foo$m1,nrow=nvar,ncol=1),
						muclus=matrix(foo$muclus,nrow=nrec+100,ncol=nvar),
						ncluster=foo$ncluster,
						psi1=matrix(foo$psi1,nrow=nvar,ncol=nvar),
						k0=foo$k0,
						sigmaclus=matrix(foo$sigmaclus,nrow=nrec+100,ncol=nuniq),
						ss=foo$ss
                      )

          randsave <- matrix(foo$randsave,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
          thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
 
          indip <- rep(0,(nvar+nvar*(nvar+1)/2+3))
 
 
          if(is.null(varnames))
          {
               varnames<-all.vars(cl)[1:nvar]
          }
 
          coeff<-NULL
          pnames1<-NULL

          if(is.null(prior$m2) && is.null(prior$s2))
		  {
	
		  }
		  else
		  {
			 for(i in 1:nvar)
			 { 
				coeff<-c(coeff,mean(thetasave[,i]))
				pnames1<-c(pnames1,paste("m1",varnames[i],sep="-"))
				indip[i]<-1
			 }
		  }
	 

         if(is.null(prior$tau1))
		 {
	
		 }
		 else
		 {
			coeff<-c(coeff,mean(thetasave[,(nvar+1)]))
			pnames1<-c(pnames1,"k0")
			indip[nvar+1]<-1
		 }
 
         if(is.null(prior$nu2))
         {
        
         }
         else
         {
			for(i in 1:(nvar*(nvar+1)/2))
			{
				coeff<-c(coeff,mean(thetasave[,(nvar+1+i)]))
				indip[(nvar+1+i)]<-1
			}
	    
             for(i in 1:nvar)
             {
                for(j in i:nvar)
                {
                  if(i==j)pnames1<-c(pnames1,paste("psi1",varnames[i],sep="-"))
                  if(i!=j)
                  {
                     tempname<-paste(varnames[i],varnames[j],sep="-")
                     pnames1<-c(pnames1,paste("psi1",tempname,sep="-"))
                  }    
                }  
             }
          }

          coeff<-c(coeff,mean(thetasave[,(nvar+nvar*(nvar+1)/2+2)]))
          pnames2<-c("ncluster")
		  indip[(nvar+nvar*(nvar+1)/2+2)]<-1

          if(alpharand==1)
          {
             coeff<-c(coeff,mean(thetasave[,(nvar+nvar*(nvar+1)/2+3)]))
             pnames2<-c(pnames2,"alpha")
	     indip[(nvar+nvar*(nvar+1)/2+3)]<-1
          }


          names(coeff)<-c(pnames1,pnames2)


          pnames1<-NULL
	  for(i in 1:nvar)
	  {
                 pnames1<-c(pnames1,paste("m1",varnames[i],sep="-"))
          }
          
          pnames1<-c(pnames1,"k0")

          for(i in 1:nvar)
          {
              for(j in i:nvar)
              {
                  if(i==j)pnames1<-c(pnames1,paste("psi1",varnames[i],sep="-"))
                  if(i!=j)
                  {
                      tempname<-paste(varnames[i],varnames[j],sep="-")
                      pnames1<-c(pnames1,paste("psi1",tempname,sep="-"))
                  }    
              }  
          }

          pnames2<-c("ncluster","alpha")
          dimnames(thetasave)<-list(NULL,c(pnames1,pnames2))

          pnamesre<-NULL              
          for(i in 1:nrec)
          {
              for(j in 1:nvar)
              {
                    tmpn<-paste("mu",varnames[j],sep="-")
                    tmpn<-paste(tmpn," (Subject=",sep="")
                    tmpn<-paste(tmpn,i,sep="")
                    tmpn<-paste(tmpn,")",sep="")
                    pnamesre<-c(pnamesre,tmpn)
              }
              
              for(j in 1:nvar)
              {
              	  for(k in j:nvar)
              	  {
                         tmpn<-paste("sigma",varnames[j],sep="-")
                         tmpn<-paste(tmpn,varnames[k],sep="-")
                         tmpn<-paste(tmpn," (Subject=",sep="")
                         tmpn<-paste(tmpn,i,sep="")
                         tmpn<-paste(tmpn,")",sep="")
                         pnamesre<-c(pnamesre,tmpn)
                     }    
              }
          }


          for(j in 1:nvar)
          {
                tmpn<-paste("mu",varnames[j],sep="-")
                tmpn<-paste(tmpn," (Prediction)",sep="")
                tmpn<-paste(tmpn,")",sep="")
                pnamesre<-c(pnamesre,tmpn)
          }
          
          for(j in 1:nvar)
          {
          	  for(k in j:nvar)
          	  {
                     tmpn<-paste("sigma",varnames[j],sep="-")
                     tmpn<-paste(tmpn,varnames[k],sep="-")
                     tmpn<-paste(tmpn," (Prediction)",sep="")
                     tmpn<-paste(tmpn,")",sep="")
                     pnamesre<-c(pnamesre,tmpn)
                 }    
          }


          for(j in 1:nvar)
          {
                tmpn<-paste(varnames[j]," (Prediction)",sep="")
                pnamesre<-c(pnamesre,tmpn)
          }

          dimnames(randsave)<-list(NULL,pnamesre)
          
          save.state <- list(thetasave=thetasave,randsave=randsave)

          x1<-NULL
          x2<-NULL
          dens<-NULL
          
          if(nvar==1)
          {
             x1 <- foo$grid1
             dens <- foo$fun1
             f <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
             grid1 <- foo$grid1
             grid2 <- foo$grid2
             fun1 <- foo$fun1
             fun2 <- foo$fun2
          }
          if(nvar==2)
          {
             x1 <- foo$grid1
             x2 <- foo$grid2
             dens <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
             f <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
             grid1 <- foo$grid1
             grid2 <- foo$grid2
             fun1 <- foo$fun1
             fun2 <- foo$fun2
          }   

          if(nvar>2)
          {
             x1 <- NULL
             x2 <- NULL
             dens <- NULL
             f <- NULL
             grid1 <- NULL
             grid2 <- NULL
             fun1 <- NULL
             fun2 <- NULL
          }   

	  z <- list(call=cl,
	            y=y,
	            varnames=varnames,
	            modelname=model.name,
	            cpo=foo$cpo,
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
                    coefficients=coeff,
                    indip=indip,
                    f=f,
                    grid1=grid1,
                    grid2=grid2,
                    fun1=fun1,
                    fun2=fun2,
                    x1=x1,
                    x2=x2,
                    dens=dens)
                 
          cat("\n\n")
 	  class(z)<-"DPdensity"
  	  return(z)
}



###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 05-10-2006.
###



"print.DPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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



"summary.DPdensity"<-function(object, hpd=TRUE, ...) 
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
    
    for(i in 1:dimen1)
    {
       if(object$indip[i]==1)
       {
           coef.p<-c(coef.p,object$coefficients[i])
           mat<-cbind(mat,thetasave[,i])
       }   
    }

    dimen1<-dim(mat)[2] 

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

    class(ans) <- "summaryDPdensity"
    return(ans)
}


"print.summaryDPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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




"plot.DPdensity"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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

   if(is(x, "DPdensity"))
   {

      if(output=="density")
      {

      # Density estimation
	
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))
	start<-(x$nrec+1)*x$nvar+(x$nrec+1)*x$nvar*(x$nvar+1)/2

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
           lines(x$grid1,x$fun1,lwd=2)
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
               lines(density(x$save.state$randsave[,(start+i)]),lwd=2)
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
                   xx<-x$save.state$randsave[,(start+i)]	    
                   yy<-x$save.state$randsave[,(start+j)]	    
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



