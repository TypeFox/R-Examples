### DPMlmm.R                   
### Fit a linear mixed model using Dirichlet Process mixture of normals
### prior for the distribution of the random effects.
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

"DPMlmm"<-
function(fixed,random,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPMlmm")


"DPMlmm.default"<-
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
         freqsub <- table(oldid)
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
         z <- matrix(1,nrow=nrec,ncol=1)
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
            startp <- 1+q
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

         tau1 <- prior$tau1
         tau2 <- prior$tau2
         tau <- c(tau1,tau2)
         if(tau1<0 || tau2<0)
         { 
                   stop("The parameters of the Gamma prior for the error variance must be possitive.\n")     
         }

         nu0 <- prior$nu0
         if(nu0<=0)
         { 
                stop("The parameter of the IW prior distribution for the normal kernel must be positive")     
         }
         tinv <- prior$tinv
         if(dim(tinv)[1]!=q || dim(tinv)[2]!=q)
         { 
                stop("Error in the dimension of the matrix of the IW prior for the covariance of the normal kernel.\n")     
         }

         if(is.null(prior$a0))
         {
            a0b0 <- c(-10,-10)
            alpha <- prior$alpha
            alphapr <- 0
         }
         else
         {
            a0b0 <- c(prior$a0,prior$b0)
            alpha <- rgamma(1,shape=prior$a0,scale=prior$b0)
            alphapr <- 1
            if(prior$a0<0 || prior$b0<0)
            { 
                   stop("The parameters of the Gamma prior for the precision parameter must be possitive.\n")     
            }
         }

         psiinv <- solve(prior$Sb)
         smu <- psiinv%*%prior$mb
         if(length(prior$mb) != q)
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

         nub <- prior$nub
         if(nub<=0)
         { 
                stop("The parameter of the IW prior for the covariance of the centering distribution must be positive")     
         }

         tbinv <- prior$tbinv
         if(dim(tbinv)[1]!=q || dim(tbinv)[2]!=q)
         { 
                stop("Error in the dimension of the matrix of the IW prior for the covariance of the centering distribution.\n")     
         }
         nu <- c(prior$nu0,prior$nub)

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
         nuniq <- q*(q+1)
         mc <- rep(0,5)         
         musave <- matrix(0,nrow=nsave,ncol=q*nsubject)
         clustsave <- matrix(0,nrow=nsave,ncol=nsubject)
         randsave <- matrix(0,nrow=nsave,ncol=q*(nsubject+1))
         thetasave <- matrix(0,nrow=nsave,ncol=q+nfixed+1+q+nuniq+2)
         cpo <- matrix(0,nrow=nrec,ncol=2)

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

				fit0 <- startlmm(fixed=resp~z-1, random = ~ z - 1 | newid, family=gaussian, q=q) 
				beta <- matrix(0,nrow=1,ncol=1)
				b <- NULL 
				for(i in 1:q)
				{
   				    b <- cbind(b,fit0$b[,i]+fit0$beta[i])
			    }
                sigma <- fit0$sigma
                sigmainv <- solve(sigma)
                sigmab <- fit0$sigma
                sigmabinv <- solve(sigmab)
				betar <- fit0$beta[1:q]
				mub <- fit0$beta[1:q]
				mu <- matrix(0,nrow=nsubject,ncol=q)
				mu[1,1:q] <- fit0$beta[1:q]
				sigma2e <- fit0$sigma2e
			}
			else
			{
				fit0 <- startlmm(fixed=resp~z+x-1, random = ~ z - 1 | newid, family=gaussian, q=q) 
				beta <- fit0$beta[(q+1):(p+q)]
				b <- NULL 
				for(i in 1:q)
				{
   				    b <- cbind(b,fit0$b[,i]+fit0$beta[i])
			    }
                sigma <- fit0$sigma
                sigmainv <- solve(sigma)
                sigmab <- fit0$sigma
                sigmabinv <- solve(sigmab)
				betar <- fit0$beta[1:q]
				mub <- fit0$beta[1:q]
				mu <- matrix(0,nrow=nsubject,ncol=q)
				mu[1,1:q] <- fit0$beta[1:q]
				sigma2e <- fit0$sigma2e 
			}
             ncluster <- 1
             ss <- rep(1,nsubject)
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
	        mu <- state$mu
	        mub <- state$mub
	        ncluster <- state$ncluster
	        sigma <- state$sigma
	        sigmab <- state$sigmab
	        sigmainv <- solve(sigma)
	        sigmabinv <- solve(sigmab)
	        sigma2e <- state$sigma2e         
	        ss <- state$ss
	        betar <- rep(0,q)
		 }
         
       #########################################################################################
       # working space
       #########################################################################################
         iflagp <- rep(0,p) 
         iflagr <- rep(0,q) 
         prob <- rep(0,nsubject+1)
         quadf <- matrix(0,nrow=q,ncol=q)
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)
         theta <- rep(0,q)
         workmhp <- rep(0,p*(p+1)/2) 
         workmp <- matrix(0,nrow=p,ncol=p) 
         workvp <- rep(0,p) 
         xty <- rep(0,p) 
         
         workmhr <- rep(0,q*(q+1)/2) 
         workmhr2 <- rep(0,q*(q+1)/2) 
         workmr <- matrix(0,nrow=q,ncol=q) 
         workvr <- rep(0,q) 
         zty <- rep(0,q) 
         ztz <- matrix(0,nrow=q,ncol=q) 
         cstrt <- matrix(0,nrow=nsubject,ncol=nsubject) 
         ccluster <- rep(0,nsubject) 

         betasave<-rep(0,(p+1))
         bsave<-matrix(0,nrow=nsubject,ncol=q)

         
       #########################################################################################
       # calling the fortran code
       #########################################################################################

         foo <- .Fortran("dpmlmm",
				maxni      =as.integer(maxni),         
				nrec       =as.integer(nrec),
				nsubject   =as.integer(nsubject),
				nfixed     =as.integer(nfixed),
				p          =as.integer(p),
				q          =as.integer(q),
				subject    =as.integer(newid),
				datastr    =as.integer(datastr), 	 	
				y          =as.double(resp),
				x          =as.double(x),	 	
				z          =as.double(z),	 
				xtx        =as.double(xtx),	 	
				a0b0       =as.double(a0b0),
				prec       =as.double(prec),	  		
				sb         =as.double(sb),	  		
				tau        =as.double(tau),	  		
				nu         =as.integer(nu),
				tinv1      =as.double(tinv),	  		 		
				smu        =as.double(smu),	  		
				psiinv     =as.double(psiinv),	  		
				tinv2      =as.double(tbinv),	  		 		
				mcmc       =as.integer(mcmcvec),
				nsave      =as.integer(nsave),
				ncluster   =as.integer(ncluster),
				ss         =as.integer(ss),
				alpha      =as.double(alpha),		
				beta       =as.double(beta),
				b          =as.double(b),			
				betar      =as.double(betar),
				mu         =as.double(mu),
				sigma2e    =as.double(sigma2e),
				sigma      =as.double(sigma),
				sigmainv   =as.double(sigmainv),
				mub        =as.double(mub),
				sigmab     =as.double(sigmab),
				sigmabinv  =as.double(sigmabinv),
				mc         =as.double(mc), 		
				cpo        =as.double(cpo),
				randsave   =as.double(randsave),
				thetasave  =as.double(thetasave),
				musave     =as.double(musave),
				clustsave  =as.integer(clustsave),
				iflagp     =as.integer(iflagp),
				workmhp    =as.double(workmhp),
				workmp     =as.double(workmp),
				workvp     =as.double(workvp),
				xty        =as.double(xty),
				iflagr     =as.integer(iflagr),
				theta      =as.double(theta),
                workmhr    =as.double(workmhr),
                workmhr2   =as.double(workmhr2),
                workmr     =as.double(workmr),
                workvr     =as.double(workvr),
				ztz        =as.double(ztz), 		
				zty        =as.double(zty), 		
				cstrt      =as.integer(cstrt),
				ccluster   =as.integer(ccluster),
				prob       =as.double(prob),
				quadf      =as.double(quadf),
				seed       =as.integer(seed),
                betasave   =as.double(betasave),
                bsave      =as.double(bsave),
				PACKAGE    ="DPpackage")	

       #########################################################################################
       # save state
       #########################################################################################

         dimen <- q+nfixed+1+q+nuniq+2
         mc <- foo$mc
         names(mc) <- c("Dbar", "Dhat", "pD", "DIC","LPML")
         thetasave <- matrix(foo$thetasave,nrow=nsave, ncol=dimen)
         randsave <- matrix(foo$randsave,nrow=nsave, ncol=q*(nsubject+1))
         musave <- matrix(foo$musave,nrow=nsave,ncol=q*nsubject)
         clustsave <- matrix(foo$clustsave,nrow=nsave,ncol=nsubject)
         
         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo <- cpom[,1]         
         fso <- cpom[,2]

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

         pnames4 <- paste("mub",colnames(z),sep="-") 	    

         pnames5 <- NULL
         for(i in 1:q)
         {
            for(j in i:q)
            {
               if(i==j) aa <-paste("sigmab",colnames(z)[i],sep="-")
               if(i!=j) aa <-paste("sigmab",colnames(z)[i],colnames(z)[j],sep="-")
               pnames5<-c(pnames5,aa)            
            }
         }

         pnames6 <- c("ncluster","alpha")

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

         colnames(thetasave)<-c(pnames1,pnames2,pnames3,pnames4,pnames5,pnames6)

         coeff<-apply(thetasave,2,mean)		

		 state <- list(alpha=foo$alpha,
	               b=matrix(foo$b,nrow=nsubject,ncol=q),
	               beta=foo$beta,
	               mu=matrix(foo$mu,nrow=nsubject,ncol=q),
	               mub=foo$mub,
	               ncluster=foo$ncluster,
	               sigma=matrix(foo$sigma,nrow=q,ncol=q),
	               sigmab=matrix(foo$sigmab,nrow=q,ncol=q),
	               sigma2e=foo$sigma2e,
	               ss=foo$ss)

		save.state <- list(thetasave=thetasave,
	                       randsave=randsave,
	                       musave=musave,
						   clustsave=clustsave)

		z <- list(	modelname=model.name,
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
					y=resp,
					possiP=possiP,
					fixed=fixed,
					mc=mc)
                 
         cat("\n\n")        

         class(z)<-c("DPMlmm")
         return(z) 
}


###                    
### Tools: anova, print, summary, plot
###
### Copyright: Alejandro Jara, 2006-2009
### Last modification: 24-04-2007.


"anova.DPMlmm"<-function(object, ...)
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



"print.DPMlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.DPMlmm"<-function(object, hpd=TRUE, ...) 
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

### Kernel variance
    dimen3<-object$nrandom*(object$nrandom+1)/2
    if(dimen3==1)
    {
       mat<-matrix(thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)],ncol=1) 
    }
    else
    {
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


    ans$kernel<-coef.table


### CPO
    ans$cpo<-object$cpo

### Baseline Information
    
    dimen4<-object$nrandom+object$nrandom*(object$nrandom+1)/2
    mat<-thetasave[,(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)]
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


    ans$base<-coef.table

### Precision parameter
    if(is.null(object$prior$a0))
    {
       dimen5<-1
       mat<-matrix(thetasave[,(dimen1+dimen2+dimen3+dimen4+1):(dimen1+dimen2+dimen3+dimen4+dimen5)],ncol=1) 
    }
    else
    {
       dimen5<-2
       mat<-thetasave[,(dimen1+dimen2+dimen3+dimen4+1):(dimen1+dimen2+dimen3+dimen4+dimen5)]
    }
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
    ans$prec<-coef.table


    coef.table<-matrix(object$mc,nrow=1,ncol=5)
    dimnames(coef.table) <- list(" ", c("Dbar", "Dhat", "pD", "DIC","LPML"))
    ans$mc<-coef.table
    
    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject

    class(ans) <- "summaryDPMlmm"
    return(ans)
}


"print.summaryDPMlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nKernel variance:\n")
    print.default(format(x$kernel, digits = digits), print.gap = 2, 
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
    else cat("No precision parameter\n")

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.DPMlmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPMlmm"))
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
			coef.p <- x$coefficients
			n <- length(coef.p)
			pnames <- names(coef.p)
			poss <- 0 
			for(i in 1:n)
			{
				if(pnames[i]==param)poss=i
			}
			if(poss==0) 
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

