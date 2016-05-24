### DPolmm.R                    
### Fit a ordinal linear mixed model with a Dirichlet Process prior
### for the random effect distribution.
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 30-04-2007.
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

"DPolmm"<-
function(fixed,random,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPolmm")


"DPolmm.default"<-
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
         if(!is.factor(resp)) stop("the response must be a factor")
         categ<-levels(resp)
         ncateg<-length(categ)
         if(ncateg <= 2) stop("the response must have 3 or more levels")
         yr <- unclass(resp)
         if(min(yr)==0)yr<-yr+1
         if(min(yr)>1)yr<-yr+min(yr)+1
         categn<-seq(1,ncateg)

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
            nburn <- 1000
            nsave <- 1000
            nskip <- 0
            ndisplay <- 100
            mcmcvec<-c(nburn,nskip,ndisplay,murand,sigmarand)
         }
         else
         {
            mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,murand,sigmarand)
            nsave<-mcmc$nsave
         }

       #########################################################################################
       # output
       #########################################################################################
         nuniq<-(q*(q+1)/2)
         mc<-rep(0,5)         
         randsave<-matrix(0,nrow=nsave,ncol=q*(nsubject+1))
         thetasave<-matrix(0,nrow=nsave,ncol=q+nfixed+q+nuniq+2+ncateg-2)
         cpo<-matrix(0,nrow=nrec,ncol=2)

       #########################################################################################
       # MLE estimation
       #########################################################################################
         
         MLEordprobit<-function(x,y,p,ncateg,nrec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta[1:p]
		cutoff<-c(0,theta[(p+1):(p+ncateg-2)])
		cc<-cumsum(c(cutoff[1], exp(cutoff[-1])))
		cutoff<-c(-100,cc,100)

		like <- pnorm(cutoff[y+1] - eta) - pnorm(cutoff[y] - eta)
                if (all(like > 0)) 
                     eval<- -sum(log(like))
                else eval<-Inf
		return(eval)
	     }
	     
	     z<-rep(0,nrec)
	     for(i in 1:nrec)
	     {
	        if(y[i]==1)
	        {
	           z[i]<-0
	        }
	        else
	        {
	           z[i]<-1
	        }
	     }
	     start<-coefficients(glm.fit(x, z, family= binomial(probit)))   
	
	     start<-c(start,seq(ncateg-2)/(ncateg-2))
	     
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }

       #########################################################################################
       # parameters depending on status
       #########################################################################################

    	 if(status==TRUE)
	 {

	        if(sigmarand==1)
	        {
                   wsigma <- prior$tinv/(prior$nu0-q-1)
                }
                else
                {
                   wsigma <- prior$sigma
                }

                bzs<-NULL
                for(i in 1:q)
                {
                   work<-rnorm(nsubject,mean=0,sd=sqrt(wsigma[i,i]))
                   bzs<-cbind(bzs,work)
                }
	 
	        if(nfixed==0)
	        {
	           beta<-matrix(0,nrow=1,ncol=1)
                   fit0<-MLEordprobit(z,yr,q,ncateg,nrec)
	           b<-matrix(0,nrow=nsubject,ncol=q)
		   bclus<-matrix(0,nrow=nsubject,ncol=q)
		   
	           for(i in 1:nsubject)
	           {
	               b[i,]<-fit0$beta[1:q]+bzs[i,]
	           }
	           bclus<-b
	           
	           if(murand==1)
	           {
  	               mu<-fit0$beta[1:q]
  	           }
  	           else
  	           {
  	               mu<-prior$mu
  	           }
  	           
	           if(sigmarand==1)
	           {
  	               sigma<-prior$tinv/(prior$nu0-q-1)
  	           }
  	           else
  	           {
  	               sigma<-prior$sigma
  	           }
	           sigmainv<-solve(sigma)
		   cutoff<-c(0,fit0$beta[(q+1):(q+ncateg-2)])
		   cutoff<-cumsum(c(cutoff[1], exp(cutoff[-1])))
	           
	        }

	        if(nfixed>0)
	        {
                   fit0<-MLEordprobit(cbind(x,z),yr,(p+q),ncateg,nrec)
	           b<-matrix(0,nrow=nsubject,ncol=q)
		   bclus<-matrix(0,nrow=nsubject,ncol=q)
                   beta<-fit0$beta[1:p]
	           for(i in 1:nsubject){
	               b[i,]<-fit0$beta[(p+1):(p+q)]+bzs[i,]
	           }
	           bclus<-b
	           
	           if(murand==1)
	           {
  	               mu<-fit0$beta[(p+1):(p+q)]
  	           }
  	           else
  	           {
  	               mu<-prior$mu
  	           }
  	           
	           if(sigmarand==1)
	           {
  	               sigma<-prior$tinv/(prior$nu0-q-1)
  	           }
  	           else
  	           {
  	               sigma<-prior$sigma
  	           }
		   sigmainv<-solve(sigma)
		   cutoff<-c(0,fit0$beta[(p+q+1):(p+q+ncateg-2)])
		   cutoff<-cumsum(c(cutoff[1], exp(cutoff[-1])))
	        }

                betar<-rep(0,q)
                ncluster<-1
                ss<-rep(1,nsubject)
	 }	
      	 if(status==FALSE)
	 {
	        alpha<-state$alpha
                b<-state$b 
                bclus<-state$bclus 
                if(nfixed>0)
                {
	           beta<-state$beta
	        }
	        else
	        {
	           beta<-rep(0,p)
	        }
	        cutoff<-state$cutoff
	        mu<-state$mu
	        ncluster<-state$ncluster
	        sigma<-state$sigma
	        sigmainv<-solve(sigma)
	        ss<-state$ss
	        betar<-rep(0,q)
	 }


       #########################################################################################
       # working space
       #########################################################################################
         ccluster<-rep(0,nsubject) 
         cstrt<-matrix(0,nrow=nsubject,ncol=nsubject)          
         iflag<-rep(0,p) 
         iflag2<-rep(0,maxni)
         iflagb<-rep(0,q) 
         prob<-rep(0,nsubject+1)
         res<-rep(0,nrec)
         seed1<-sample(1:29000,1)
         seed2<-sample(1:29000,1)
         seed<-c(seed1,seed2)
         theta<-rep(0,q)
         work1<-matrix(0,nrow=p,ncol=p)
         workb1<-matrix(0,nrow=q,ncol=q)
         workb2<-matrix(0,nrow=q,ncol=q)
         workmh1<-rep(0,p*(p+1)/2) 
         workmh2<-rep(0,q*(q+1)/2) 
         workmh3<-rep(0,q*(q+1)/2) 
         workk1<-matrix(0,nrow=maxni,ncol=q) 
         workkv1<-rep(0,maxni) 
         workkm1<-matrix(0,nrow=maxni,ncol=maxni) 
         workkm2<-matrix(0,nrow=maxni,ncol=maxni) 
         workv1<-rep(0,p) 
         workvb1<-rep(0,q) 
         xty<-rep(0,p) 
         y<-rep(0,nrec) 
         ywork<-rep(0,maxni) 
         zty<-rep(0,q) 
         ztz<-matrix(0,nrow=q,ncol=q) 

         betasave<-rep(0,(p+ncateg-1))
         bsave<-matrix(0,nrow=nsubject,ncol=q)
         
       #########################################################################################
       # calling the fortran code
       #########################################################################################
         foo <- .Fortran("dpolmmp",
         	datastr    =as.integer(datastr),
 	 	maxni      =as.integer(maxni),         
 	 	ncateg     =as.integer(ncateg),         
 	 	nrec       =as.integer(nrec),
 	 	nsubject   =as.integer(nsubject),
 	 	nfixed     =as.integer(nfixed),
 	 	p          =as.integer(p),
 	 	q          =as.integer(q),
 	 	subject    =as.integer(newid),
 		x          =as.double(x),	 	
 		xtx        =as.double(xtx),	 	
 		y          =as.double(y),
 		yr         =as.integer(yr),
 		z          =as.double(z),	 
 		a0b0       =as.double(a0b0),
 		nu0        =as.integer(nu0),
 		prec       =as.double(prec),	 
 		psiinv     =as.double(psiinv),	  		
 		sb         =as.double(sb),	  		
 		smu        =as.double(smu),	  		
 		tinv       =as.double(tinv),	  		 		
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
 		randsave   =as.double(randsave),
 		thetasave  =as.double(thetasave),
 		cpo        =as.double(cpo),
 		alpha      =as.double(alpha),		
 		b          =as.double(b),		
                bclus      =as.double(bclus),		
 		beta       =as.double(beta),
 		betar      =as.double(betar),
 		cutoff     =as.double(cutoff),
 		mu         =as.double(mu),
 		ncluster   =as.integer(ncluster),
 		sigma      =as.double(sigma),
 		ss         =as.integer(ss),
 		mc         =as.double(mc), 		
 		ccluster   =as.integer(ccluster),
 		cstrt      =as.integer(cstrt),
 		iflag      =as.integer(iflag),
 		iflag2     =as.integer(iflag2),
 		iflagb     =as.integer(iflagb),
 		prob       =as.double(prob),
 		res        =as.double(res),
 		seed       =as.integer(seed),
 		sigmainv   =as.double(sigmainv),
 		theta      =as.double(theta),
 		work1      =as.double(work1),
 		workb1     =as.double(workb1),
 		workb2     =as.double(workb2),
                workmh1    =as.double(workmh1),
                workmh2    =as.double(workmh2),
                workmh3    =as.double(workmh3),
                workk1     =as.double(workk1),
 		workkv1    =as.double(workkv1),
 		workkm1    =as.double(workkm1),
 		workkm2    =as.double(workkm2),
 		workv1     =as.double(workv1),
 		workvb1    =as.double(workvb1),
 		xty        =as.double(xty),
 		ywork      =as.double(ywork),
 		zty        =as.double(zty), 		
 		ztz        =as.double(ztz), 
                betasave   =as.double(betasave),
                bsave      =as.double(bsave),
		PACKAGE    ="DPpackage")	

       #########################################################################################
       # save state
       #########################################################################################

         dimen<-q+nfixed+q+nuniq+2+ncateg-2
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
 	 
         pnames2 <- paste("mu",colnames(z),sep="-") 	     	 
         pnames3<-NULL
         for(i in 1:q)
         {
            for(j in i:q)
            {
               if(i==j) aa <-paste("sigma",colnames(z)[i],sep="-")
               if(i!=j) aa <-paste("sigma",colnames(z)[i],colnames(z)[j],sep="-")
               pnames3<-c(pnames3,aa)            
            }
         }
         
         pnames4<-paste("cutoff",2:(ncateg-1))
         
         pnames5<- c("ncluster","alpha")
         
         colnames(thetasave)<-c(pnames1,pnames2,pnames3,pnames4,pnames5)
         
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
        
	 model.name<-"Bayesian semiparametric ordinal linear mixed effect model"		

         coeff<-apply(thetasave, 2, mean)

	 state <- list(alpha=foo$alpha,
	               b=matrix(foo$b,nrow=nsubject,ncol=q),
	               bclus=matrix(foo$bclus,nrow=nsubject,ncol=q),
	               beta=foo$beta,
	               cutoff=foo$cutoff,
	               mu=foo$mu,
	               ncluster=foo$ncluster,
	               sigma=matrix(foo$sigma,nrow=q,ncol=q),
	               ss=foo$ss)

	 save.state <- list(thetasave=thetasave,randsave=randsave)

	 z<-list(modelname=model.name,
	         coefficients=coeff,call=cl,
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
                 ncateg=ncateg,
                 nsave=nsave,
                 y=resp,
                 possiP=possiP,
                 fixed=fixed,
                 mc=mc)
                 
         cat("\n\n")        

         class(z)<-c("DPolmm")
         return(z) 
}


###
### Tools for DPolmm: anova, print, summary, plot
###
### Copyright: Alejandro Jara, 2006-2007
### Last modification: 01-04-2007.

"anova.DPolmm"<-function(object, ...)
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



"print.DPolmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"summary.DPolmm"<-function(object, hpd=TRUE, ...) 
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
    
    ans <- c(object[c("call", "modelname")])

    ans$coefficients<-coef.table


### CPO
    ans$cpo<-object$cpo


### Baseline Information

    dimen2<-object$nrandom+object$nrandom*(object$nrandom+1)/2
    
    if(dimen2==1)
    {
       mat<-matrix(thetasave[,(dimen1+1):(dimen1+dimen2)],ncol=1) 
    }
    else
    {
       mat<-thetasave[,(dimen1+1):(dimen1+dimen2)]
    }
    
    coef.p<-object$coefficients[(dimen1+1):(dimen1+dimen2)]
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


### Cutoff Points

    dimen3<-object$ncateg-2
    
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

    ans$cutoff<-coef.table


### Precision parameter

    if(is.null(object$prior$a0))
    {
      dimen4<-1
      coef.p<-object$coefficients[(dimen1+dimen2+dimen3+1)]
      mat<-matrix(thetasave[,(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)],ncol=1)
    }
    else
    {
      dimen4<-length(object$coefficients)-(dimen1+dimen2+dimen3)
      coef.p<-object$coefficients[(dimen1+dimen2+dimen3+1):length(object$coefficients)]
      mat<-thetasave[,(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)]

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

    coef.table<-matrix(object$mc,nrow=1,ncol=5)
    dimnames(coef.table) <- list(" ", c("Dbar", "Dhat", "pD", "DIC","LPML"))
    ans$mc<-coef.table

    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject

    class(ans) <- "summaryDPolmm"
    return(ans)
}


"print.summaryDPolmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nThreshold:\n")
    print.default(format(x$cutoff, digits = digits), print.gap = 2, 
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


"plot.DPolmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPolmm"))
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
           
           if(is.null(x$prior$a0))
           {
               cat("")
           }
           else
           {
               title1<-paste("Trace of",pnames[n],sep=" ")
               title2<-paste("Density of",pnames[n],sep=" ")       
               plot(ts(x$save.state$thetasave[,n]),main=title1,xlab="MCMC scan",ylab=" ")
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
   }

}


