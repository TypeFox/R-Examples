### PTlmm.R                   
### Fit a linear mixed model using a Mxture of Multivariate Polya Trees 
### prior for the distribution of the random effects.
###
### Copyright: Alejandro Jara and Tim Hanson, 2007-2012.
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
###      Tim Hanson
###      Division of Biostatistics
###      University of Minnesota
###      School of Public Health
###      A460 Mayo Building, 
###      MMC 303
###      420 Delaware St SE
###      Minneapolis, MN 55455
###      Voice: 612-626-7075  URL  : http://www.biostat.umn.edu/~hanson/
###      Fax  : 612-626-0660  Email: hanson@biostat.umn.edu
###

"PTlmm"<-
function(fixed,random,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PTlmm")


"PTlmm.default"<-
function(fixed,
         random,
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
         keep <- is.element(nm, c("data", "na.action"))
         for (i in nm[!keep]) m[[i]] <- NULL
         
         allvars <- c(all.vars(fixed), all.vars(random))

         Terms <- if (missing(data)) 
              terms(fixed)
         else terms(fixed, data = data)
         
         cl$fixed <- eval(fixed)
         cl$random <- eval(random)
         m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
         environment(m$formula) <- environment(fixed)
         m$drop.unused.levels <- TRUE
         m[[1]] <- as.name("model.frame")
         mf <- eval.parent(m)

       #########################################################################################
       # data structure
       #########################################################################################
         nrec <- dim(mf)[1]
         resp <- mf[,1]

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
         xtx<-t(x)%*%x
         
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
            prec1<-matrix(0,nrow=1,ncol=1)
            sb<-matrix(0,nrow=1,ncol=1)
            b0<-rep(0,1)
         }
         else
         {
            b0<-prior$beta0
            prec1<-solve(prior$Sbeta0)
            sb<-prec1%*%b0

            if(length(b0)!=p)
            { 
                   stop("Error in the dimension of the mean of the normal prior for the fixed effects.\n")     
            }

            if(dim(prec1)[1]!=p || dim(prec1)[2]!=p)
            { 
                   stop("Error in the dimension of the covariance of the normal prior for the fixed effects.\n")     
            }

         }

         tau1<-prior$tau1
         tau2<-prior$tau2
         tau<-c(tau1,tau2)
         if(tau1<0 || tau2<0)
         { 
                   stop("The parameters of the Gamma prior for the error variance must be possitive.\n")     
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
			murand <- 0 
            if(is.null(prior$mu))
            { 
               stop("The vector *mu* must be specified in the prior object when it is not considered as random.\n")     
            }
            if(length(prior$mu) != q)
            { 
               stop("Error in the dimension of the mean the centering distribution.\n")     
            }
            mu <- prior$mu
            prec2<-diag(1,q)
            mu0<-rep(0,q)
          }
          else
          {
            murand <- 1
			prec2 <- solve(prior$Sb)
			mu0 <- prec2%*%prior$mub

            if(length(prior$mub) != q)
            { 
               stop("Error in the dimension of the mean of the normal prior for the mean of the centering distribution.\n")     
            }

            if(is.null(dim(prec2)) && q==1) 
            {
               stop("Error in the dimension of the covariance of the normal prior for the mean of the centering distribution.\n")     
            }   

            if(!is.null(dim(prec2)) && ( dim(prec2)[1]!=q || dim(prec2)[2]!=q ))
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
            sigma <- prior$sigma
            nu0 <- -1
            tinv<-diag(1,q)
            a0b0<-c(a0b0,nu0,tau1,tau2)
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
			a0b0<-c(a0b0,nu0,tau1,tau2)
		 }
         
		 if(is.null(prior$frstlprob))
		 {
            frstlprob<-0
		 }
		 else
		 {
            frstlprob<-0
            if(prior$frstlprob)frstlprob<-1
		 }

         if(is.null(prior$typepr))
         {
            typepr <- 1
         }
         else
         {
            typepr <- prior$typepr
         }
         if(q==1)
         {
            typepr <- 0
            ortho <- diag(q)
         }

       #########################################################################################
       # mcmc specification
       #########################################################################################
         if(missing(mcmc))
         {
            tune1 <- -1.1
            tune2 <- 1.1
            tune3 <- 1.1
            tune4 <- -1.1
            M <- log(nsubject,2^q)
            nburn <- 1000
            nsave <- 1000
            nskip <- 0
            ndisplay <- 100
            nbase <- 1
            samplef <- 1 
            mcmcvec <- c(nburn,nskip,ndisplay,nbase,M)
            a0b0 <- c(a0b0,tune1,tune2,tune3,tune4)
         }
         else
         {
            if(is.null(mcmc$tune1))
            {  
               tune1 <- -1.1
            }
            else
            {
               tune1 <- mcmc$tune1
            }

            if(is.null(mcmc$tune2))
            {
               tune2 <- -1.1
            }
            else
            {
               tune2 <- mcmc$tune2
            }

            if(is.null(mcmc$tune3))
            {
               tune3 <- -1.1
            }
            else
            {
               tune3 <- mcmc$tune3
            }

            if(is.null(mcmc$tune4))
            {
               tune4 <- -1.1
            }
            else
            {
               tune4 <- mcmc$tune4
            }


            if(is.null(mcmc$nbase))
            {
               nbase <- 1
            }
            else
            {
               nbase <- mcmc$nbase
            }

            if(is.null(prior$M))
            {
               M <- log(nsubject,q)
            } 
            else
            {
               M <- prior$M
            }

            if(is.null(mcmc$samplef))
            {
               samplef <- 1
            } 
            else
            {
               samplef <- mcmc$samplef
            }

            mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,nbase,M)
            nsave <- mcmc$nsave
            a0b0 <- c(a0b0,tune1,tune2,tune3,tune4)
         }

       #########################################################################################
       # output
       #########################################################################################
         acrate <- rep(0,5)
         nuniq <- q*(q+1)
         cpo <- matrix(0,nrow=nrec,ncol=2)
         mc <- rep(0,5)
         randsave <- matrix(0,nrow=nsave,ncol=q*(nsubject+1))
         thetasave <- matrix(0,nrow=nsave,ncol=q+nfixed+1+q+nuniq+1+q*q)


       #########################################################################################
       # parameters depending on status
       #########################################################################################
		 startlmm <- function(fixed,random,family,q)
         {
			 library(nlme)
			 library(MASS)
			 fit0 <- glmmPQL(fixed=fixed, random=random, family=family, verbose = FALSE) 
			 beta <- fit0$coeff$fixed
			 b <- fit0$coeff$random$newid
			 sigma <- getVarCov(fit0)[1:q,1:q]
			 sigma2e <- fit0$sigma^2
			 out <- list(beta=beta,b=b,sigma=sigma,sigma2e=sigma2e)
			 return(out)
         }

    	 if(status==TRUE)
   	     {
            if(nfixed==0)
            {
				fit0 <- startlmm(fixed=resp~z-1, random = ~ z - 1 | newid, family=gaussian,q=q) 
				beta <- matrix(0,nrow=1,ncol=1)
				b <- NULL 
				for(i in 1:q)
				{
   				    b <- cbind(b,fit0$b[,i]+fit0$beta[i])
			    }
			}
			else
			{

				fit0 <- startlmm(fixed=resp~z+x-1, random = ~ z - 1 | newid, family=gaussian,q=q) 
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
			betar<-rep(0,q)
			ortho <- matrix(rnorm(q*q),nrow=q,ncol=q)
			sigma2e <- fit0$sigma2e

		 }	
      	 if(status==FALSE)
		 {
	        alpha <- state$alpha
			b <- state$b 
			if(nfixed>0)
			{
	           beta <- state$beta
	        }
	        else
	        {
	           beta <- rep(0,p)
	        }
			if(murand==1)
			{
   	            mu <- state$mu
   	        }    
			if(sigmarand==1)
			{
   	            sigma <- state$sigma
   	        }    
	        sigma2e <- state$sigma2e    
			if(typepr==1)
			{
				ortho <- state$ortho 
			}
			else
			{
				ortho <- matrix(rnorm(q*q),nrow=q,ncol=q)
			}                 
	        betar <-rep(0,q)	        
		 }
         
       #########################################################################################
       # working space
       #########################################################################################
         narea <- 2^q
         curr <- c(alpha,sigma2e)
         bz <- matrix(0,nrow=nsubject,ncol=q)
         bzc <- matrix(0,nrow=nsubject,ncol=q)
         iflagp <- rep(0,p) 
         iflagr <- rep(0,q) 
         limw <- rep(0,q)
         linf <- rep(0,q)
         lsup <- rep(0,q)
         massi <- rep(0,narea)
         parti <- rep(0,q)
         pattern <- rep(0,q)
         propvr <- matrix(0,nrow=q,ncol=q)
         res <- rep(0,nrec)
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         sigmac <- matrix(0,nrow=q,ncol=q)
         sigmainv <- matrix(0,nrow=q,ncol=q)
         sigmainvc <- matrix(0,nrow=q,ncol=q)
         theta <- rep(0,q)
         thetac <- rep(0,q)
         whicho <- rep(0,nsubject)
         whichn <- rep(0,nsubject)      
         workmhp1 <- rep(0,p*(p+1)/2) 
         workmhr <- rep(0,q*(q+1)/2) 
         workmhr2 <- rep(0,q*(q+1)/2) 
         workmp1 <- matrix(0,nrow=p,ncol=p)
         workmr <- matrix(0,nrow=q,ncol=q) 
         workmr1 <- matrix(0,nrow=q,ncol=q) 
         workmr2 <- matrix(0,nrow=q,ncol=q) 
         workvp1 <- rep(0,p) 
         workvr <- rep(0,q)
         xty <- rep(0,p) 
         ybar <- rep(0,narea)
         
         betasave <- rep(0,(p+1))
         bsave <- matrix(0,nrow=nsubject,ncol=q)

         mcmcvec <- c(mcmcvec,seed1,seed2,typepr,murand,sigmarand,frstlprob,samplef)

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         foo <- .Fortran("ptlmm",
				datastr    =as.integer(datastr),
				maxni      =as.integer(maxni),
				nrec       =as.integer(nrec),
				nsubject   =as.integer(nsubject),
				nfixed     =as.integer(nfixed),
				p          =as.integer(p),
				q          =as.integer(q),
				subject    =as.integer(newid),
				x          =as.double(x),
				xtx        =as.double(xtx),	 	
				y          =as.double(resp),
				z          =as.double(z),
				a0b0       =as.double(a0b0),
				mu0        =as.integer(mu0),
				prec1      =as.double(prec1),	 
				prec2      =as.double(prec2),	 
				sb         =as.double(sb),	  		
				tinv       =as.double(tinv),	  		 		
				mcmc       =as.integer(mcmcvec),
				nsave      =as.integer(nsave),
                acrate     =as.double(acrate),   
				cpo        =as.double(cpo),
				randsave   =as.double(randsave),
				thetasave  =as.double(thetasave),
				curr       =as.double(curr),
				b          =as.double(b),		
				beta       =as.double(beta),
				betar      =as.double(betar),
				mu         =as.double(mu),
				sigma      =as.double(sigma),
                ortho      =as.double(ortho),
				mc         =as.double(mc),
				iflagp     =as.integer(iflagp),
				res        =as.double(res),
				workmp1    =as.double(workmp1),
                workmhp1   =as.double(workmhp1),
				workvp1    =as.double(workvp1),
				xty        =as.double(xty),
				iflagr     =as.integer(iflagr),
				parti      =as.integer(parti),
				whicho     =as.integer(whicho),
				whichn     =as.integer(whichn),
				bz         =as.double(bz),
				bzc        =as.double(bzc),
				limw       =as.double(limw),
				linf       =as.double(linf),
				lsup       =as.double(lsup),
				propvr     =as.double(propvr),
				sigmainv   =as.double(sigmainv),
				theta      =as.double(theta),
				thetac     =as.double(thetac),
				workmhr    =as.double(workmhr),
				workmr     =as.double(workmr),
				workmr1    =as.double(workmr1),
				workmr2    =as.double(workmr2),
				workvr     =as.double(workvr),
				ybar       =as.double(ybar),
				sigmac     =as.double(sigmac),
				sigmainvc  =as.double(sigmainvc),
				workmhr2   =as.double(workmhr2),
				massi      =as.integer(massi),
				pattern    =as.integer(pattern),
                betasave   =as.double(betasave),
                bsave      =as.double(bsave),
				PACKAGE    ="DPpackage")	
		

       #########################################################################################
       # save state
       #########################################################################################

         dimen<-q+nfixed+1+q+nuniq+1+q*q
         mc<-foo$mc
         names(mc)<-c("Dbar", "Dhat", "pD", "DIC","LPML")
         thetasave<-matrix(foo$thetasave,nrow=nsave, ncol=dimen)
         randsave<-matrix(foo$randsave,nrow=nsave, ncol=q*(nsubject+1))

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
 	 
		 pnames2<-"residual"
 	 
         pnames3 <- paste("mu",colnames(z),sep="-") 	    

         pnames4 <- NULL
         for(i in 1:q)
         {
            for(j in i:q)
            {
               if(i==j) aa <-paste("sigma",colnames(z)[i],sep="-")
               if(i!=j) aa <-paste("sigma",colnames(z)[i],colnames(z)[j],sep="-")
               pnames4<-c(pnames4,aa)            
            }
         }
         
         pnames5<- c("alpha",paste("ortho",seq(q*q)))

         pnames6 <- NULL
         for(i in 1:q)
         {
            for(j in i:q)
            {
               if(i==j) aa <-paste("R.E.Cov",colnames(z)[i],sep="-")
               if(i!=j) aa <-paste("R.E.Cov",colnames(z)[i],colnames(z)[j],sep="-")
               pnames6<-c(pnames6,aa)            
            }
         }
        
         colnames(thetasave)<-c(pnames1,pnames2,pnames3,pnames4,pnames5,pnames6)
         
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
         
         
		 model.name<-"Bayesian semiparametric linear mixed effect model"		

         coeff<-apply(thetasave,2,mean)
		 names(coeff)<-c(pnames1,pnames2,pnames3,pnames4,pnames5,pnames6)

		 state <- list(	alpha=foo$curr[1],
						b=matrix(foo$b,nrow=nsubject,ncol=q),
						beta=foo$beta,
						mu=foo$mu,
						sigma=matrix(foo$sigma,nrow=q,ncol=q),
						sigma2e=foo$curr[2],
						ortho=matrix(foo$ortho,nrow=q,ncol=q))

		 save.state <- list(	thetasave=thetasave,
							randsave=randsave)

         acrate<-foo$acrate

		z<-list(	modelname=model.name,
					coefficients=coeff,
					call=cl,
					prior=prior,
					mcmc=mcmc,
					state=state,
					save.state=save.state,
					nrec=foo$nrec,
					nsubject=foo$nsubject,
					nfixed=foo$nfixed,
					nrandom=foo$q,
					cpo=cpo,
					fso=fso,
					alphapr=alphapr,
					namesre1=namesre,
					namesre2=colnames(z),
					z=z,
					x=x,
					mf=mf,
					dimen=dimen,
					acrate=acrate,
					possiP=possiP,
					fixed=fixed,
					m=M,
					murand=murand,
					sigmarand=sigmarand,
					frstlprob=frstlprob,
					mc=mc,
					typepr=typepr,
					samplef=samplef)
                 
         cat("\n\n")        

         class(z)<-c("PTlmm")
         return(z) 
}


###                    
### Tools: anova, print, summary, plot
###
### Copyright: Alejandro Jara, 2006-2007-2008
### Last modification: 30-04-2007.


"anova.PTlmm"<-function(object, ...)
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



"print.PTlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.PTlmm"<-function(object, hpd=TRUE, ...) 
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

### Fixed part of the model

    dimen1<-object$nrandom+object$nfixed
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

### Residual variance
    
    dimen2<-1
    mat<-matrix(thetasave[,(dimen1+dimen2)],ncol=1)     
    coef.p<-object$coefficients[(dimen1+dimen2)]
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

    ans$resvar<-coef.table


### CPO
    ans$cpo<-object$cpo

### Baseline Information
    
    dimen3<-object$nrandom+object$nrandom*(object$nrandom+1)/2
    mat<-thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]
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
    ans$base<-coef.table

### Precision parameter
    dimen4<-1

    if(!is.null(object$prior$a0))
    {
       mat<-matrix(thetasave[,(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)],ncol=1) 
       coef.p<-object$coefficients[(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)]
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
    }

### Partition
    dimen5<-object$nrandom*object$nrandom

    if(object$typepr==1)
    {
       mat<-matrix(thetasave[ ,(dimen1+dimen2+dimen3+dimen4+1):(dimen1+dimen2+dimen3+dimen4+dimen5)],ncol=dimen5) 
       coef.p<-object$coefficients[(dimen1+dimen2+dimen3+dimen4+1):(dimen1+dimen2+dimen3+dimen4+dimen5)]
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
    
       ans$part<-coef.table
    }

### Random effects covariance
    if(object$samplef==1)
    {
    dimen6<-object$nrandom*(object$nrandom+1)/2
    if(dimen6==1)
    {
       mat<-matrix(thetasave[,(dimen1+dimen2+dimen3+dimen4+dimen5+1):(dimen1+dimen2+dimen3+dimen4+dimen5+dimen6)],ncol=1) 
    }
    else
    {
       mat<-thetasave[,(dimen1+dimen2+dimen3+dimen4+dimen5+1):(dimen1+dimen2+dimen3+dimen4+dimen5+dimen6)]
    }
    coef.p<-object$coefficients[(dimen1+dimen2+dimen3+dimen4+dimen5+1):(dimen1+dimen2+dimen3+dimen4+dimen5+dimen6)]
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
    ans$recov<-coef.table
    }


### Model comparison
    coef.table<-matrix(object$mc,nrow=1,ncol=5)
    dimnames(coef.table) <- list(" ", c("Dbar", "Dhat", "pD", "DIC","LPML"))
    ans$mc<-coef.table
    
    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject

    ans$acrate<-object$acrate

    class(ans) <- "summaryPTlmm"
    return(ans)
}


"print.summaryPTlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nResidual variance:\n")
    print.default(format(x$resvar, digits = digits), print.gap = 2, 
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
    else cat("\nNo precision parameter\n")

    if (length(x$part)) {
        cat("\nPartition:\n")
        print.default(format(x$part, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$recov)) {
    cat("\nRandom effects variance:\n")
    print.default(format(x$recov, digits = digits), print.gap = 2, 
          quote = FALSE)
    }
    else cat("\nFunctional parameters were not sampled\n")

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.PTlmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "PTlmm"))
   {
        if(is.null(param))
        {
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:n)
           {
               if(is.null(x$prior$a0) && pnames[i]=="alpha")
               {
               
               }
               if((x$typepr==0) && pnames[i]=="typep")
               {
               
               }
               else
               {
                  title1<-paste("Trace of",pnames[i],sep=" ")
                  title2<-paste("Density of",pnames[i],sep=" ")       
                  plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
                  fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }   
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
            fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
        }
   }

}

