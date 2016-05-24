### DPglmm.R                   
### Fit a generalized linear mixed model with a Dirichlet Process prior
### for the random effect distribution
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 25-09-2009.
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

"DPglmm"<-
function(fixed,random,family,offset,n,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPglmm")


"DPglmm.default"<-
function(fixed,
         random,
         family,
         offset=NULL,
         n=NULL,
         prior,
         mcmc,
         state,
         status, 
         data=sys.frame(sys.parent()),
         na.action=na.fail)
{
       #########################################################################################
       # call parameters
       #########################################################################################
         m <- mcall <- cl <- match.call()
         nm <- names(m)[-1]
         keep <- is.element(nm, c("data", "na.action","offset"))
         for (i in nm[!keep]) m[[i]] <- NULL
         
         if(is.null(n))
         {  
            allvars <- c(all.vars(fixed), all.vars(random))
         }
         else
         {
            allvars <- c(all.vars(fixed), all.vars(random), "n")         
         }   
         
         Terms <- if (missing(data)) 
              terms(fixed)
         else terms(fixed, data = data)
         
         off <- attr(Terms, "offset")
         if (length(off <- attr(Terms, "offset"))) 
            allvars <- c(allvars, as.character(attr(Terms, "variables"))[off + 1])

         cl$fixed <- eval(fixed)
         cl$random <- eval(random)
         m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
         environment(m$formula) <- environment(fixed)
         m$drop.unused.levels <- TRUE
         m[[1]] <- as.name("model.frame")
         mf <- eval.parent(m)

       #########################################################################################
       # checking input
       #########################################################################################
         if(is.null(n))
         {
              collapsed <- 0         
         }
         else
         {
              collapsed <- 1
              if(family$family!="binomial" || family$link!="logit")
              { 
                   stop("The total number of cases only can be used in the logit link .\n")     
              }
         }

       #########################################################################################
       # data structure
       #########################################################################################
         nrec <- dim(mf)[1]
         resp <- mf[,1]
         roffset <- model.offset(mf)
         if (is.null(roffset)) roffset <-rep(0,nrec)

         crandom <- all.vars(random)
         namesre <- (names(mf)==crandom[length(crandom)])
         oldid <- mf[,namesre]
         freqsub<-table(oldid)
         namesre <- names(freqsub)         
         nsubject <- length(namesre)
         newid <- seq(1,nrec)
         for(i in 1:nsubject)
         {
             newid[oldid==namesre[i]] <- i
         }
         
         if(is.null(n))
         {
           ntrials <- rep(1,nrec)
         }
         else
         {
           ntrials <- (names(mf)=="n")
           ntrials <- mf[,ntrials]
         }  

         maxni <- max(freqsub)
         idrec <- seq(1,nrec)
         datastr <- matrix(0,nrow=nsubject,ncol=maxni+1)
         datastr[,1] <- freqsub
         for(i in 1:nsubject)
         {
             for(j in 1:freqsub[i])
             {
                 datastr[i,(j+1)] <- idrec[newid==i][j] 
             }
         }

       #########################################################################################
       # model structure
       #########################################################################################
         q <- length(crandom)
         z<-matrix(1,nrow=nrec,ncol=1)
         colnames(z) <- "(Intercept)"
         nvarrand <- "(Intercept)"

         if(q>1)
         {
            for(i in 1:(q-1))
            {
               zwork <- (names(mf)==crandom[i])
               zwork <- matrix(mf[,zwork],nrow=nrec,ncol=1)
               colnames(zwork) <- crandom[i]
               nvarrand <- c(nvarrand,crandom[i])
               z <- cbind(z,zwork)
            }
         }   
         
         cfixed <- all.vars(fixed)
         cfixed <- cfixed[-1]

         for(i in 1:q)
         {
            if(sum(nvarrand[i]==cfixed) != 0)
            {
               stop("Covariates cannot be included in the random and fixed part of the model")
            }   
         }   
         x <- model.matrix(fixed,data=mf)
         p <- dim(x)[2]
         x <- x[,-1]
         p <- p-1
         nfixed <- p
         
         if(p==0)
         {
            nfixed <- 0
            p <- 1
            x <- matrix(0,nrow=nrec,ncol=1)
         }
         
       #########################################################################################
       # elements for Pseudo Countour Probabilities' computation
       #########################################################################################
         possiP <- NULL
         if(nfixed>0)
         {
            mat <- attr(Terms,"factors")
            namfact <- colnames(mat)
            nvar <- dim(mat)[1]
            nfact <- dim(mat)[2]
            possiP <- matrix(0,ncol=2,nrow=nfact)
            if (missing(data)) dataF <- model.frame(formula=fixed,xlev=NULL)
               dataF <- model.frame(formula=fixed,data,xlev=NULL)
            namD <- names(dataF)
            isF <- sapply(dataF, function(x) is.factor(x) || is.logical(x))
            nlevel <- rep(0,nvar)
            for(i in 1:nvar)
            {
                if(isF[i])
                {
                   nlevel[i]<-length(table(dataF[[i]]))
                }
                else
                {
                   nlevel[i]<-1
                }
            }
            startp<-1+q
            for(i in 1:nfact)
            {
                tmp1<-1
                for(j in 1:nvar)
                {
                    if(mat[j,i]==1 && isF[j])
                    {
                       tmp1<-tmp1*(nlevel[j]-1)
                    }
                }
                endp<-startp+tmp1-1
                possiP[i,1]<-startp    
                possiP[i,2]<-endp
                startp<-endp+1
            }
            dimnames(possiP)<-list(namfact,c("Start","End"))
         }   

       #########################################################################################
       # prior information
       #########################################################################################
         if(nfixed==0)
         {
            prec<-matrix(0,nrow=1,ncol=1)
            sb<-matrix(0,nrow=1,ncol=1)
            b0<-rep(0,1)
         }
         else
         {
            b0<-prior$beta0
            prec<-solve(prior$Sbeta0)
            sb<-prec%*%b0

            if(length(b0)!=p)
            { 
                   stop("Error in the dimension of the mean of the normal prior for the fixed effects.\n")     
            }

            if(dim(prec)[1]!=p || dim(prec)[2]!=p)
            { 
                   stop("Error in the dimension of the covariance of the normal prior for the fixed effects.\n")     
            }

         }

         if(family$family=="Gamma")
         {
            tau1<-prior$tau1
            tau2<-prior$tau2
            tau<-c(tau1,tau2)
            if(tau1<0 || tau2<0)
            { 
                   stop("The parameters of the Gamma prior for the dispersion parameter must be possitive.\n")     
            }
         }    

         if(is.null(prior$a0))
         {
            a0b0<-c(-10,-10)
            alpha<-prior$alpha
            alphapr<-0
         }
         else
         {
            a0b0<-c(prior$a0,prior$b0)
            alpha<-rgamma(1,shape=prior$a0,scale=prior$b0)
            alphapr<-1
            if(prior$a0<0 || prior$b0<0)
            { 
                   stop("The parameters of the Gamma prior for the precision parameter must be possitive.\n")     
            }
         }

		if(is.null(prior$mub))
		{
			murand<-0 
            if(is.null(prior$mu))
            { 
               stop("The vector *mu* must be specified in the prior object when it is not considered as random.\n")     
            }
            if(length(prior$mu) != q)
            { 
               stop("Error in the dimension of the mean the centering distribution.\n")     
            }
            psiinv<-diag(1,q)
            smu<-rep(0,q)
         }
         else
         {
            murand<-1
			psiinv<-solve(prior$Sb)
			smu<-psiinv%*%prior$mub
            if(length(prior$mub) != q)
            { 
               stop("Error in the dimension of the mean of the normal prior for the mean of the centering distribution.\n")     
            }

            if(is.null(dim(psiinv)) && q==1) 
            {
               stop("Error in the dimension of the covariance of the normal prior for the mean of the centering distribution.\n")     
            }   

            if(!is.null(dim(psiinv)) && ( dim(psiinv)[1]!=q || dim(psiinv)[2]!=q ))
            { 
               stop("Error in the dimension of the covariance of the normal prior for the mean of the centering distribution.\n")     
            }
         }

		if(is.null(prior$nu0))
		{
            sigmarand<-0
            if(is.null(prior$sigma))
            { 
               stop("The matrix *sigma* must be specified in the prior object when it is not considered as random.\n")     
            }
            nu0 <- -1
            tinv<-diag(1,q)
		}
		else
		{
			sigmarand<-1
			nu0<-prior$nu0
			if(nu0<=0)
			{ 
                stop("The parameter of the IW prior distribution must be positive")     
			}

			tinv<-prior$tinv
			if(dim(tinv)[1]!=q || dim(tinv)[2]!=q)
			{ 
                stop("Error in the dimension of the matrix of the IW prior for the covariance of the centering distribution.\n")     
			}
         }

       #########################################################################################
       # mcmc specification
       #########################################################################################
         if(missing(mcmc))
         {
            tune4 <- 1.1
            nburn <- 1000
            nsave <- 1000
            nskip <- 0
            ndisplay <- 100
            mpar <- 1
            mcmcvec<-c(nburn,nskip,ndisplay,murand,sigmarand)
         }
         else
         {
            if(is.null(mcmc$tune1))
            {
               tune4 <- 1.1
            }
            else
            {
               tune4 <- mcmc$tune1
            }
            mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,murand,sigmarand)
            nsave<-mcmc$nsave

            if(is.null(mcmc$mpar))
            {
               mpar <- 1
            }
            else
            {
               mpar <- mcmc$mpar
            }   
         }


       #########################################################################################
       # output
       #########################################################################################
         nuniq <- (q*(q+1)/2)
         acrate <- rep(0,2)
         cpo <- matrix(0,nrow=nrec,ncol=2)
         dispp <- 0
         if(family$family=="Gamma") dispp <- 1
         mc <- rep(0,5)                           
         randsave <- matrix(0,nrow=nsave,ncol=q*(nsubject+1))
         thetasave <- matrix(0,nrow=nsave,ncol=q+nfixed+dispp+q+nuniq+2)

       #########################################################################################
       # parameters depending on status
       #########################################################################################
		 startglmm <- function(fixed,random,family,q)
         {
			 library(nlme)
			 library(MASS)
			 fit0 <- glmmPQL(fixed=fixed, random=random, family=family, verbose = FALSE) 
			 beta <- fit0$coeff$fixed
			 b <- fit0$coeff$random$newid
			 sigma <- getVarCov(fit0)[1:q,1:q]
			 out <- list(beta=beta,b=b,sigma=sigma)
			 return(out)
         }

         if(status==TRUE)
         {
            resp2 <- resp
            if(family$family=="binomial")
            {
                if(family$link=="logit")
                {
                    if(!is.null(n)) resp2<-cbind(resp,ntrials-resp)
                }
            }   
			 
            if(nfixed==0)
            {
				fit0 <- startglmm(fixed=resp2~z-1+offset(roffset), random = ~ z - 1 | newid, family=family,q=q) 
				beta <- matrix(0,nrow=1,ncol=1)
				b <- NULL 
				for(i in 1:q)
				{
   				    b <- cbind(b,fit0$b[,i]+fit0$beta[i])
			    }
			}
			else
			{

				fit0 <- startglmm(fixed=resp2~z+x-1+offset(roffset), random = ~ z - 1 | newid, family=family,q=q) 
				beta <- fit0$beta[(q+1):(p+q)]
				b <- NULL 
				for(i in 1:q)
				{
   				    b <- cbind(b,fit0$b[,i]+fit0$beta[i])
			    }
			}

			if(sigmarand==1)
			{
				sigma <- fit0$sigma
			}
			else
			{
				sigma <- prior$sigma
			}
		 
		    if(murand==1)
			{
				mu <- fit0$beta[1:q]
			}
			else
			{
				mu <- prior$mu
			}
			bclus <- b 
			betar <- rep(0,q)
			ncluster <- nsubject
			ss <- seq(1,nsubject)
			if(family$family=="Gamma") tau <- c(tau,1.1,tune4)
			sigmainv <- solve(sigma)
          }	
         if(status==FALSE)
         {
             alpha <- state$alpha
             b <- state$b 
             bclus <- state$bclus 
             if(nfixed>0)
             {
				 beta <- state$beta
			 }
			 else
			 {
				 beta <- rep(0,p)
			 }
             mu <- state$mu
             ncluster <- state$ncluster
             sigma <- state$sigma
             ss <- state$ss
             if(family$family=="Gamma")tau <- c(tau,1/state$phi,tune4)

             betar <- rep(0,q)
             sigmainv <- solve(sigma)
         }

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         if(family$family=="binomial")
         {
            if(family$link=="logit")
            {
                   resp <- cbind(resp,ntrials)
                   betac <- rep(0,p)
                   cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)
                   ccluster <- rep(0,nsubject) 
                   iflag <- rep(0,p)
                   iflagb <- rep(0,q) 
                   prob <- rep(0,(nsubject+1))
                   quadf <- matrix(0,nrow=q,ncol=q)
                   seed1 <- sample(1:29000,1)
                   seed2 <- sample(1:29000,1)
                   seed <- c(seed1,seed2)
                   theta <- rep(0,q)
                   thetac <- rep(0,q)
                   workb1 <- matrix(0,nrow=q,ncol=q)
                   workb2 <- matrix(0,nrow=q,ncol=q)
                   workmh1 <- rep(0,p*(p+1)/2) 
                   workmh2 <- rep(0,q*(q+1)/2) 
                   workmh3 <- rep(0,q*(q+1)/2) 
                   workv1 <- rep(0,p) 
                   workvb1 <- rep(0,q) 
                   workvb2 <- rep(0,q) 
                   xtx <- matrix(0,nrow=p,ncol=p) 
                   xty <- rep(0,p) 
                   zty <- rep(0,q) 
                   ztz <- matrix(0,nrow=q,ncol=q) 

                   betasave <- rep(0,p)
                   bsave <- matrix(0,nrow=nsubject,ncol=q)

                 # fitting the model
                 
                   foo <- .Fortran("dpglmmlogita",
                        datastr   =as.integer(datastr),
                        maxni     =as.integer(maxni),         
                        nrec      =as.integer(nrec),
                        nsubject  =as.integer(nsubject),
                        nfixed    =as.integer(nfixed),
                        p         =as.integer(p),
                        q         =as.integer(q),
                        subject   =as.integer(newid),
                        x         =as.double(x),	 	
                        y         =as.integer(resp),
                        z         =as.double(z),	 
                        a0b0      =as.double(a0b0),
                        b0        =as.double(b0),
                        nu0       =as.integer(nu0),
                        prec      =as.double(prec),	 
                        psiinv    =as.double(psiinv),	  		
                        sb        =as.double(sb),	  		
                        smu       =as.double(smu),	  		
                        tinv      =as.double(tinv),	  		 		
                        mcmc      =as.integer(mcmcvec),
                        nsave     =as.integer(nsave),
                        acrate    =as.double(acrate),
                        cpo       =as.double(cpo),
                        randsave  =as.double(randsave),
                        thetasave =as.double(thetasave),
                        alpha     =as.double(alpha),		
                        b         =as.double(b),		
                        bclus     =as.double(bclus),		
                        beta      =as.double(beta),
                        betar     =as.double(betar),
                        mu        =as.double(mu),
                        ncluster  =as.integer(ncluster),
                        sigma     =as.double(sigma),
                        sigmainv  =as.double(sigmainv),
                        ss        =as.integer(ss),
						mc        =as.double(mc),                        
                        betac     =as.double(betac),
                        cstrt     =as.integer(cstrt),
                        ccluster  =as.integer(ccluster),
                        iflag     =as.integer(iflag),
                        iflagb    =as.integer(iflagb),
                        prob      =as.double(prob),
                        quadf     =as.double(quadf),
                        seed      =as.integer(seed),
                        theta     =as.double(theta),
                        thetac    =as.double(thetac),
                        workb1    =as.double(workb1),
                        workb2    =as.double(workb2),
                        workmh1   =as.double(workmh1),
                        workmh2   =as.double(workmh2),
                        workmh3   =as.double(workmh3),
                        workv1    =as.double(workv1),
                        workvb1   =as.double(workvb1),
                        workvb2   =as.double(workvb2),
                        xtx       =as.double(xtx),
                        xty       =as.double(xty),
                        zty       =as.double(zty), 		
                        ztz       =as.double(ztz), 		
                        betasave   =as.double(betasave),
                        bsave      =as.double(bsave),
                        PACKAGE   ="DPpackage")	            
            }

            if(family$link=="probit")
            {
                # specifying the working space
                  xtx <- t(x)%*%x
                  lat <- rep(0,nrec)
                  ccluster <- rep(0,nsubject) 
                  cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)                            
                  iflag <- rep(0,p) 
                  iflag2 <- rep(0,maxni)
                  iflagb <- rep(0,q) 
                  prob <- rep(0,(nsubject+1))
                  quadf <- matrix(0,nrow=q,ncol=q)
                  res <- rep(0,nrec)
                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)
                  theta <- rep(0,q)
                  work1 <- matrix(0,nrow=p,ncol=p)
                  workb1 <- matrix(0,nrow=q,ncol=q)
                  workb2 <- matrix(0,nrow=q,ncol=q)
                  workmh1 <- rep(0,p*(p+1)/2) 
                  workmh2 <- rep(0,q*(q+1)/2) 
                  workmh3 <- rep(0,q*(q+1)/2) 
                  workk1 <- matrix(0,nrow=maxni,ncol=q) 
                  workkv1 <- rep(0,maxni) 
                  workkm1 <- matrix(0,nrow=maxni,ncol=maxni) 
                  workkm2 <- matrix(0,nrow=maxni,ncol=maxni) 
                  workv1 <- rep(0,p) 
                  workvb1 <- rep(0,q) 
                  workvb2 <- rep(0,q) 
                  xty <- rep(0,p) 
                  ywork <- rep(0,maxni) 
                  zty <- rep(0,q) 
                  ztz <- matrix(0,nrow=q,ncol=q) 

                  betasave <- rep(0,p)
                  bsave <- matrix(0,nrow=nsubject,ncol=q)
                
                # fitting the model
                
                  foo <- .Fortran("dpglmmprob",
                        datastr   =as.integer(datastr),
                        maxni     =as.integer(maxni),         
                        nrec      =as.integer(nrec),
                        nsubject  =as.integer(nsubject),
                        nfixed    =as.integer(nfixed),
                        p         =as.integer(p),
                        q         =as.integer(q),
                        subject   =as.integer(newid),
                        x         =as.double(x),	 	
                        xtx       =as.double(xtx),
                        y         =as.double(lat), 
                        yr        =as.integer(resp),
                        z         =as.double(z),	 
                        a0b0      =as.double(a0b0),
                        nu0       =as.integer(nu0),
                        prec      =as.double(prec),	 
                        psiinv    =as.double(psiinv),	  		
                        sb        =as.double(sb),	  		
                        smu       =as.double(smu),	  		
                        tinv      =as.double(tinv),	  		 		
                        mcmc      =as.integer(mcmcvec),
                        nsave     =as.integer(nsave),
                        randsave  =as.double(randsave),
                        thetasave =as.double(thetasave),
                        cpo       =as.double(cpo),		
                        alpha     =as.double(alpha),		
                        b         =as.double(b),		
                        bclus     =as.double(bclus),		
                        beta      =as.double(beta),
                        betar     =as.double(betar),
                        mu        =as.double(mu),
                        ncluster  =as.integer(ncluster),
                        sigma     =as.double(sigma),
                        ss        =as.integer(ss),
						mc        =as.double(mc),                        
                        ccluster  =as.integer(ccluster),
						cstrt     =as.integer(cstrt), 		                        
                        iflag     =as.integer(iflag),
                        iflag2    =as.integer(iflag2),
                        iflagb    =as.integer(iflagb),
                        prob      =as.double(prob),
                        quadf     =as.double(quadf),
                        res       =as.double(res),
                        seed      =as.integer(seed),
                        sigmainv  =as.double(sigmainv),
                        theta     =as.double(theta),
                        work1     =as.double(work1),
                        workb1    =as.double(workb1),
                        workb2    =as.double(workb2),
                        workmh1   =as.double(workmh1),
                        workmh2   =as.double(workmh2),
                        workmh3   =as.double(workmh3),
                        workk1    =as.double(workk1),
                        workkv1   =as.double(workkv1),
                        workkm1   =as.double(workkm1),
                        workkm2   =as.double(workkm2),
                        workv1    =as.double(workv1),
                        workvb1   =as.double(workvb1),
                        workvb2   =as.double(workvb2),
                        xty       =as.double(xty),
                        ywork     =as.double(ywork),
                        zty       =as.double(zty), 		
                        ztz       =as.double(ztz), 		
                        betasave  =as.double(betasave),
                        bsave     =as.double(bsave),
                        PACKAGE   ="DPpackage")	            
            }
         }

         if(family$family=="poisson")
         {
            if(family$link=="log")
            {
  	        # specifying the working space
                  betac <- rep(0,p)
                  ccluster <- rep(0,nsubject) 
                  cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)
                  iflagp <- rep(0,p)
                  iflagb <- rep(0,q) 
                  prob <- rep(0,(nsubject+mpar))
                  newtheta <- matrix(0,nrow=q,ncol=mpar)
                  quadf <- matrix(0,nrow=q,ncol=q)
                  seed <- c(sample(1:29000,1),sample(1:29000,1))
                  theta <- rep(0,q)
                  thetac <- rep(0,q)
                  workp1 <- matrix(0,nrow=p,ncol=p)
                  workp2 <- matrix(0,nrow=p,ncol=p)
                  workb1 <- matrix(0,nrow=q,ncol=q)
                  workb2 <- matrix(0,nrow=q,ncol=q)
                  workmh1 <- rep(0,(p*(p+1)/2)) 
                  workmh2 <- rep(0,(q*(q+1)/2)) 
                  workmh3 <- rep(0,(q*(q+1)/2)) 
                  workvp1 <- rep(0,p) 
                  workvb1 <- rep(0,q) 
                  workvb2 <- rep(0,q) 
                  xtx <- matrix(0,nrow=p,ncol=p) 
                  xty <- rep(0,p) 
                  zty <- rep(0,q) 
                  ztz <- matrix(0,nrow=q,ncol=q) 

                  betasave <- rep(0,p)
                  bsave <- matrix(0,nrow=nsubject,ncol=q)
                
                  foo <- .Fortran("dpglmmpois",
                        datastr   =as.integer(datastr),
                        maxni     =as.integer(maxni),         
                        mpar      =as.integer(mpar),
                        nrec      =as.integer(nrec),
                        nsubject  =as.integer(nsubject),
                        nfixed    =as.integer(nfixed),
                        p         =as.integer(p),
                        q         =as.integer(q),
                        subject   =as.integer(newid),
                        x         =as.double(x),	 	
                        y         =as.integer(resp),
                        z         =as.double(z),
                        roffset   =as.double(roffset),
                        a0b0      =as.double(a0b0),
                        b0        =as.double(b0),
                        nu0       =as.integer(nu0),
                        prec      =as.double(prec),	 
                        psiinv    =as.double(psiinv),	  		
                        sb        =as.double(sb),	  		
                        smu       =as.double(smu),	  		
                        tinv      =as.double(tinv),	  		 		
                        mcmc      =as.integer(mcmcvec),
                        nsave     =as.integer(nsave),
                        acrate    =as.double(acrate),
                        cpo       =as.double(cpo),
                        randsave  =as.double(randsave),
                        thetasave =as.double(thetasave),
                        alpha     =as.double(alpha),		
                        b         =as.double(b),		
                        bclus     =as.double(bclus),		
                        beta      =as.double(beta),
                        betar     =as.double(betar),
                        mu        =as.double(mu),
                        ncluster  =as.integer(ncluster),
                        sigma     =as.double(sigma),
                        sigmainv  =as.double(sigmainv),
                        ss        =as.integer(ss),
						mc        =as.double(mc),                        
                        betac     =as.double(betac),
                        ccluster  =as.integer(ccluster),
                        iflagp    =as.integer(iflagp),
                        iflagb    =as.integer(iflagb),
                        newtheta  =as.double(newtheta),
                        prob      =as.double(prob),
                        quadf     =as.double(quadf),
                        seed      =as.integer(seed),
                        theta     =as.double(theta),
                        thetac    =as.double(thetac),
                        workp1    =as.double(workp1),
                        workp2    =as.double(workp2),
                        workb1    =as.double(workb1),
                        workb2    =as.double(workb2),
                        workmh1   =as.double(workmh1),
                        workmh2   =as.double(workmh2),
                        workmh3   =as.double(workmh3),
                        workvp1   =as.double(workvp1),
                        workvb1   =as.double(workvb1),
                        workvb2   =as.double(workvb2),
                        xtx       =as.double(xtx),
                        xty       =as.double(xty),
                        zty       =as.double(zty), 		
                        ztz       =as.double(ztz), 		
                        cstrt     =as.integer(cstrt),
                        betasave  =as.double(betasave),
                        bsave     =as.double(bsave),
                        PACKAGE   ="DPpackage")	  
            }
         }   

         if(family$family=="Gamma")
         {
            if(family$link=="log")
            {
  	        # specifying the working space
                  acrate <- rep(0,3)
                  betac <- rep(0,p)
                  ccluster <- rep(0,nsubject) 
                  cstrt <- matrix(0,nrow=nsubject,ncol=nsubject)                  
                  iflag <- rep(0,p)
                  iflagb <- rep(0,q) 
                  prob <- rep(0,nsubject+2)
                  quadf <- matrix(0,nrow=q,ncol=q)
                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)
                  theta <- rep(0,q)
                  thetac <- rep(0,q)
                  workb1 <- matrix(0,nrow=q,ncol=q)
                  workb2 <- matrix(0,nrow=q,ncol=q)
                  workmh1 <- rep(0,p*(p+1)/2) 
                  workmh2 <- rep(0,q*(q+1)/2) 
                  workmh3 <- rep(0,q*(q+1)/2) 
                  workv1 <- rep(0,p) 
                  workvb1 <- rep(0,q) 
                  workvb2 <- rep(0,q) 
                  xtx <- matrix(0,nrow=p,ncol=p) 
                  xty <- rep(0,p) 
                  zty <- rep(0,q) 
                  ztz <- matrix(0,nrow=q,ncol=q) 

                  betasave<-rep(0,(p+1))
                  bsave<-matrix(0,nrow=nsubject,ncol=q)

                # fitting the model

                  foo <- .Fortran("dpglmmgam",
                        datastr   =as.integer(datastr),
                        maxni     =as.integer(maxni),         
                        nrec      =as.integer(nrec),
                        nsubject  =as.integer(nsubject),
                        nfixed    =as.integer(nfixed),
                        p         =as.integer(p),
                        q         =as.integer(q),
                        subject   =as.integer(newid),
                        x         =as.double(x),	 	
                        y         =as.double(resp),
                        z         =as.double(z),
                        roffset   =as.double(roffset),
                        a0b0      =as.double(a0b0),
                        b0        =as.double(b0),
                        nu0       =as.integer(nu0),
                        prec      =as.double(prec),	 
                        psiinv    =as.double(psiinv),	  		
                        sb        =as.double(sb),	  		
                        smu       =as.double(smu),
                        tau       =as.double(tau),
                        tinv      =as.double(tinv),	  		 		
                        mcmc      =as.integer(mcmcvec),
                        nsave     =as.integer(nsave),
                        acrate    =as.double(acrate),
                        cpo       =as.double(cpo),
                        randsave  =as.double(randsave),
                        thetasave =as.double(thetasave),
                        alpha     =as.double(alpha),		
                        b         =as.double(b),		
                        bclus     =as.double(bclus),		
                        beta      =as.double(beta),
                        betar     =as.double(betar),
                        mu        =as.double(mu),
                        ncluster  =as.integer(ncluster),
                        sigma     =as.double(sigma),
                        sigmainv  =as.double(sigmainv),
                        ss        =as.integer(ss),
						mc        =as.double(mc),
                        betac     =as.double(betac),
                        cstrt     =as.integer(cstrt),
                        ccluster  =as.integer(ccluster),
                        iflag     =as.integer(iflag),
                        iflagb    =as.integer(iflagb),
                        prob      =as.double(prob),
                        quadf     =as.double(quadf),
                        seed      =as.integer(seed),
                        theta     =as.double(theta),
                        thetac    =as.double(thetac),
                        workb1    =as.double(workb1),
                        workb2    =as.double(workb2),
                        workmh1   =as.double(workmh1),
                        workmh2   =as.double(workmh2),
                        workmh3   =as.double(workmh3),
                        workv1    =as.double(workv1),
                        workvb1   =as.double(workvb1),
                        workvb2   =as.double(workvb2),
                        xtx       =as.double(xtx),
                        xty       =as.double(xty),
                        zty       =as.double(zty), 		
                        ztz       =as.double(ztz), 		
                        betasave  =as.double(betasave),
                        bsave     =as.double(bsave),
                        PACKAGE   ="DPpackage")	  
            }
         }   

       #########################################################################################
       # save state
       #########################################################################################

         mc<-foo$mc
         names(mc)<-c("Dbar", "Dhat", "pD", "DIC","LPML")

         dimen <- q+nfixed+dispp+q+nuniq+2
         thetasave <- matrix(foo$thetasave,nrow=nsave, ncol=dimen)
         randsave <- matrix(foo$randsave,nrow=nsave, ncol=q*(nsubject+1))

         cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo<-cpom[,1]         
         fso<-cpom[,2]
 
         if(nfixed==0)
         {
            pnames1 <- colnames(z)
         }
         if(nfixed==1)
         {
            pnames1 <- c(colnames(z),cfixed)
         }
         if(nfixed>1)
         {
            pnames1 <- c(colnames(z),colnames(x))
         }   

         pnames2 <- paste("mu",colnames(z),sep="-") 	    

         pnames3 <- NULL
         for(i in 1:q)
         {
            for(j in i:q)
            {
               if(i==j) aa <-paste("sigma",colnames(z)[i],sep="-")
               if(i!=j) aa <-paste("sigma",colnames(z)[i],colnames(z)[j],sep="-")
               pnames3<-c(pnames3,aa)            
            }
         }

         pnames4 <- c("ncluster","alpha")

         if(family$family=="Gamma")
         {
            colnames(thetasave) <- c(pnames1,"phi",pnames2,pnames3,pnames4)
         }
         else
         {
            colnames(thetasave) <- c(pnames1,pnames2,pnames3,pnames4)
         }   

         qnames <- NULL
         for(i in 1:nsubject)
         {
             for(j in 1:q)
             {
                 idname <- paste("(Subject",namesre[i],sep="=")
                 idname <- paste(idname,")")
                 qnamestemp <- paste(colnames(z)[j],idname,sep=" ")
                 qnames <- c(qnames,qnamestemp)
             }
         }
         for(j in 1:q)
         {
             qnamestemp <- paste("pred",colnames(z)[j],sep="-")
             qnames <- c(qnames,qnamestemp)
         }
         colnames(randsave) <- qnames
         
         model.name <- "Bayesian semiparametric generalized linear mixed effect model"		


         coeff<-apply(thetasave,2,mean)		

         state <- list(alpha=foo$alpha,
	               b=matrix(foo$b,nrow=nsubject,ncol=q),
	               bclus=matrix(foo$bclus,nrow=nsubject,ncol=q),
	               beta=foo$beta,
	               mu=foo$mu,
	               ncluster=foo$ncluster,
	               sigma=matrix(foo$sigma,nrow=q,ncol=q),
	               ss=foo$ss,
	               phi=1/foo$tau[3])

         save.state <- list(thetasave=thetasave,randsave=randsave)

         z<-list(modelname=model.name,
				 coefficients=coeff,
				 call=cl,
                 prior=prior,
                 mcmc=mcmc,
                 state=state,
                 save.state=save.state,
                 nrec=foo$nrec,
                 nsubject=
                 foo$nsubject,
                 nfixed=foo$nfixed,
                 nrandom=foo$q,
                 cpo=cpo,
                 fso=fso,
                 alphapr=alphapr,
                 prior=prior,
                 namesre1=namesre,
                 namesre2=colnames(z),
                 z=z,
                 x=x,
                 mf=mf,
                 dimen=dimen,
                 acrate=foo$acrate,
                 y=resp,
                 dispp=dispp,
                 possiP=possiP,
                 fixed=fixed,
                 mc=mc)
                 
         cat("\n\n")        

         class(z)<-c("DPglmm")
         return(z) 
}


###                    
### Tools: anova, print, summary, plot
###
### Copyright: Alejandro Jara, 2006-2007
### Last modification: 25-03-2007.

"anova.DPglmm"<-function(object, ...)
{

######################################################################################
cregion<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a simultaneous credible region for a vector 
#  parameter from the MCMC sample
# 
#  Reference: Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995)
#             Bayesian computation and stochastic systems (with Discussion)
#             Statistical Science, vol. 10, 3 - 66, page 30
#  and        Held, L. (2004) Simultaneous inference in risk assessment; a Bayesian 
#             perspective In: COMPSTAT 2004, Proceedings in Computational 
#             Statistics (J. Antoch, Ed.) 213 - 222, page 214
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are computed.
######################################################################################
{
    #Basic information
     nmonte<-dim(x)[1]
     p<-dim(x)[2]
     
    #Ranks for each component
     ranks <- apply(x, 2, rank, ties.method="first")
     
    #Compute the set S={max(nmonte+1-min r_i(t) , max r_i(t)): t=1,..,nmonte}
     left <- nmonte + 1 - apply(ranks, 1, min)
     right <- apply(ranks, 1, max)
     S <- apply(cbind(left, right), 1, max)
     S <- S[order(S)]
    
    #Compute the credible region
     k <- floor(nmonte*probs)     
     tstar <- S[k]
     out<-list()
     for(i in 1:length(tstar))
     {
        upelim <- x[ranks == tstar[i]]
        lowlim <- x[ranks == nmonte + 1 - tstar[i]]    
        out[[i]] <- rbind(lowlim, upelim)
        rownames(out[[i]]) <- c("Lower", "Upper")
        colnames(out[[i]]) <- colnames(x)
     }
     names(out) <- paste(probs)
     return(out)
}

######################################################################################
cint<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a credible interval from the MCMC sample
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are to be computed.
######################################################################################
{
    #Compute the credible interval
     delta<-(1-probs)/2
     lprobs<-cbind(delta,probs+delta) 
     out<-matrix(quantile(x,probs=lprobs),ncol=2)
     colnames(out) <- c("Lower","Upper")
     rownames(out) <- paste(probs)
     return(out)
}

######################################################################################
hnulleval<-function(mat,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     npar<-dim(mat)[2]   
     lower<-rep(0,npar)
     upper<-rep(0,npar)
     for(i in 1:npar)
     {
        lower[i]<-mat[1,i]< hnull[i]
        upper[i]<-mat[2,i]> hnull[i]
     }
     total<-lower+upper
     out<-(sum(total==2) == npar)
     return(out)
}

######################################################################################
hnulleval2<-function(vec,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     lower<-vec[1]< hnull
     upper<-vec[2]> hnull

     total<-lower+upper
     out<-(total==2)
     return(out)
}


######################################################################################
pcp<-function(x,hnull=NULL,precision=0.001,prob=0.95,digits=digits)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Region)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-rep(0,dim(x)[2])
    if (dim(x)[2]!=length(hnull)) stop("Dimension of x and hnull must be equal!!")

    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cregion(x,probs=probsf)

    is.hnull <- hnulleval(cr[[2]],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval(cr[[length(cr)]],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval(cr[[i]],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1 - probs[first]
       }
    }
    output <- list(cr=cr[[1]], prob=prob, pval=pval,hnull=hnull)
    return(output)
}


######################################################################################
pcp2<-function(x,hnull=NULL,precision=0.001,prob=0.95)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Interval)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-0
    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cint(x,probs=probsf)

    is.hnull <- hnulleval2(cr[2,],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval2(cr[(neval+1),],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval2(cr[i,],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1-probs[first]
       }
    }
    output <- list(cr=cr[1,], prob=prob, pval=pval,hnull=hnull)
    return(output)
}

######################################################################################
######################################################################################
######################################################################################
    if(object$nfixed>0)
    {
       possiP<-object$possiP
       nfact<-dim(possiP)[1]
       P<-rep(0,nfact)
       df<-rep(0,nfact)
    
       for(i in 1:nfact)
       {
           df[i]<-1
           if((possiP[i,2]-possiP[i,1])>0)
           { 
              x<-matrix(object$save.state$thetasave[,possiP[i,1]:possiP[i,2]])
              foo<-pcp(x=x) 
              P[i]<-foo$pval
              df[i]<-(possiP[i,2]-possiP[i,1])+1
           }
           else
           {
              x<-object$save.state$thetasave[,possiP[i,1]:possiP[i,2]]
              foo<-pcp2(x=x) 
              P[i]<-foo$pval
           }
       }

       table <- data.frame(df,P) 
       dimnames(table) <- list(rownames(possiP), c("Df","PsCP"))
       structure(table, heading = c("Table of Pseudo Contour Probabilities\n", 
        paste("Response:", deparse(formula(object$fixed)[[2]]))), class = c("anovaPsCP",
        "data.frame"))
    }    
}



"print.DPglmm" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) {
        cat("Posterior Inference of Parameters:\n")
        if(x$alphapr==1){
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)}
        if(x$alphapr==0){
        print.default(format(x$coefficients[1:(length(x$coefficients)-1)], digits = digits), print.gap = 2, 
            quote = FALSE)}

    }
    else cat("No coefficients\n")

    if(!is.null(x$acrate))
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   
   
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.DPglmm"<-function(object, hpd=TRUE, ...) 
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

### Fixed part of the model

    dimen1<-object$nrandom+object$nfixed+object$dispp
    if(dimen1==1)
    {
       mat<-matrix(thetasave[,1:dimen1],ncol=1) 
    }
    else
    {
       mat<-thetasave[,1:dimen1]
    }

    coef.p<-object$coefficients[1:dimen1]
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

    ans <- c(object[c("call", "modelname")])
    ans$coefficients<-coef.table


### CPO
    ans$cpo<-object$cpo


### Baseline Information
    
    dimen2<-object$nrandom+object$nrandom*(object$nrandom+1)/2
    mat<-thetasave[,(dimen1+1):(dimen1+dimen2)]
    coef.p<-object$coefficients[(dimen1+1):(dimen1+dimen2)]
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

    ans$base<-coef.table


### Precision parameter
    if(is.null(object$prior$a0))
    {
       dimen3<-1
       mat<-matrix(thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)],ncol=1) 
    }
    else
    {
       dimen3<-2
       mat<-thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]
    }
    coef.p<-object$coefficients[(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]
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
    ans$prec<-coef.table

    coef.table<-matrix(object$mc,nrow=1,ncol=5)
    dimnames(coef.table) <- list(" ", c("Dbar", "Dhat", "pD", "DIC","LPML"))
    ans$mc<-coef.table
    
    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject
    ans$acrate<-object$acrate

    class(ans) <- "summaryDPglmm"
    return(ans)
}


"print.summaryDPglmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 

    cat("\nModel's performance:\n")
    print.default(format(x$mc, digits = digits), print.gap = 2, 
    quote = FALSE)
            
    if (length(x$coefficients)) {
        cat("\nRegression coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

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
    else cat("No precision parameter\n")


    if(!is.null(x$acrate))
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}



"plot.DPglmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPglmm"))
   {
        if(is.null(param))
        {
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)
           
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
            coef.p<-x$coefficients
	    n<-length(coef.p)
	    pnames<-names(coef.p)
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
            if(param=="ncluster")
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


