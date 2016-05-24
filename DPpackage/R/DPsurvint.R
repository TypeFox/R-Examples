### DPsurvint.R                    
### Fit a semiparametric aft model for interval censored data.
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 22-12-2006.
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

"DPsurvint"<-
function(formula,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPsurvint")

"DPsurvint.default"<-
function(formula,
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

	 cl <- match.call()
	 mf <- match.call(expand.dots = FALSE)
	 m <- match(c("formula", "data","na.action"), names(mf), 0)
	 mf <- mf[c(1, m)]
	 mf$drop.unused.levels <- TRUE
	 mf[[1]] <- as.name("model.frame")
	 mf <- eval(mf, parent.frame())

         #########################################################################################
         # data structure
         #########################################################################################

  	 y<- model.response(mf,"numeric")
  	 nrec<-length(y[,1])
  	 x<-model.matrix(formula)
  	 namesxm<-colnames(x)
  	 p<-dim(x)[2]
  	 x<-as.matrix(x[,2:p])
  	 p<-(p-1)

         type<-rep(2,nrec)
         for(i in 1:nrec){
	    type[y[,1]==-999]<-1
	    type[y[,2]==-999]<-3
         }

         #########################################################################################
         # Elements for Pseudo Countour Probabilities' computation
         #########################################################################################

         tt<-terms(formula,data=data)
         mat<-attr(tt,"factors")
         namfact<-colnames(mat)
         nvar<-dim(mat)[1]
         nfact<-dim(mat)[2]
         possiP<-matrix(0,ncol=2,nrow=nfact)
         dataF<-model.frame(formula,data,xlev=NULL)
         namD<-names(dataF)
         isF<-sapply(dataF, function(x) is.factor(x) || is.logical(x))
         nlevel<-rep(0,nvar)
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
         startp<-1
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


         #########################################################################################
         # prior information
         #########################################################################################

	 if(is.null(prior$a0))
	 {
		a0b0<-c(-1,-1)
		alpha<-prior$alpha
	 }
	 else
	 {
	 	a0b0<-c(prior$a0,prior$b0)
	 	alpha<-rgamma(1,prior$a0,prior$b0)
	 }
	
	 betapm<-prior$beta0
	 betapv<-prior$Sbeta0
	 m0<-prior$m0
	 s0<-prior$s0
	 tau<-c(prior$tau1,prior$tau2)
	
         #########################################################################################
         # mcmc specification
         #########################################################################################

         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave
         
  	 tmpy<- y[type==2,]
         mid<-matrix((tmpy[,1]+tmpy[,2])/2,ncol=1)
         nmid<-dim(mid)[1]

         tmpx<- x[type==2,]
         tmpx<- cbind(rep(1,nmid),tmpx)
         beta<-solve(t(tmpx)%*%tmpx)%*%t(tmpx)%*%log(mid)

         s2<-sum((log(mid)-tmpx%*%beta)**2)/(nmid-p-1)
	 propv<-solve(t(tmpx)%*%tmpx)*s2
	 
	 beta<--beta[2:(p+1)]
	 propv<-propv[2:(p+1),2:(p+1)]
	 
         if(is.null(mcmc$tune))
         {
            tune=1
         }
         else
         {
            tune<-mcmc$tune
         }

	 propv<-diag(tune,p)%*%solve(solve(betapv)+solve(propv))%*%diag(tune,p)

         #########################################################################################
         # output
         #########################################################################################

	 thetasave <- matrix(0, nrow=mcmc$nsave, ncol=p+4)
	 randsave  <- matrix(0, nrow=mcmc$nsave, ncol=nrec+1)
	 
	
         #########################################################################################
         # parameters depending on status
         #########################################################################################

	 if(status==TRUE)
	 {
	 	eta <- x %*% beta
		v <- rep(0,nrec)
		
		for(i in 1:nrec)
		{
		   if(type[i]==1)v[i]<- y[i,2]*exp(eta[i])/2
		   if(type[i]==2)v[i]<-(y[i,1]*exp(eta[i])+y[i,2]*exp(eta[i]))/2
		   if(type[i]==3)v[i]<-y[i,1]*exp(eta[i])+1
		}
		
		mu<-0
		sigma<-s2
	 }	
	 if(status==FALSE)
	 {
	        alpha<-state$alpha 
		beta<-state$beta
		v<-state$v
		eta <- x %*% beta
		mu<-state$mu
		sigma<-state$sigma
	 }

         #########################################################################################
         # working space
         #########################################################################################

         extra<-1 
	 maxint<-4*nrec+1
	 maxend<-4*nrec
	 ntotal<-nrec
	 maxm<-as.integer(log(0.0001)/log((ntotal+1)/(ntotal+2)))

 	 acrate<-0
	 betac<-rep(0,p)
	 clusts<-matrix(0,nrow=nrec,ncol=(nrec+1))	 
	 cpo<-rep(0,nrec)
         etan<-rep(0,nrec)
	 endp<-rep(0,maxend)
	 endp2<-rep(0,maxend)
 	 iflag<-rep(0,p)
 	 imaxs<-rep(0,nrec)
 	 imaxsc<-rep(0,nrec)
 	 imins<-rep(0,nrec)
 	 iminsc<-rep(0,nrec)
 	 index<-rep(0,maxm)	 
	 intcount<-rep(0,maxint)
	 intcount2<-rep(0,maxint)
	 intind<-rep(0,(nrec+1))
	 intind2<-rep(0,(nrec+1))
	 intposso<-matrix(0,nrow=maxint,ncol=(nrec+1))
         intpossn<-matrix(0,nrow=maxint,ncol=(nrec+1))
         limr<-matrix(0,nrow=nrec,ncol=2)         
         linfs<-rep(0,nrec)
         linfsc<-rep(0,nrec)
         lsups<-rep(0,nrec)
         lsupsc<-rep(0,nrec)
	 mass<-rep(0,maxint)
	 ncluster<-0
	 prob<-rep(0,maxint)
	 s<-rep(0,nrec)	 
	 seed<-c(sample(1:29000,1),sample(1:29000,1),sample(1:29000,1))
	 uvec<-rep(0,maxm)
	 vvec<-rep(0,maxm)
	 vnew<-rep(0,(nrec+1))
	 vnew2<-rep(0,nrec)
         wvec<-rep(0,maxm)
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)


         #########################################################################################
         # calling the fortran code
         #########################################################################################
 
 	 foo <- .Fortran("dpsurvint",
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		x         =as.double(x),	 	
		y         =as.double(y),
		interind  =as.integer(type),
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(mcmc$nsave),
		propv     =as.double(propv),
		extra     =as.integer(extra),		
		a0b0      =as.integer(a0b0),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		m0        =as.double(m0),		
		s0        =as.double(s0),		
		tau       =as.double(tau),
		acrate    =as.double(acrate),
		thetasave =as.double(thetasave),
		randsave  =as.double(randsave),
                cpo       =as.double(cpo),
                ncluster  =as.integer(ncluster),
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		mu        =as.double(mu),
		sigma2    =as.double(sigma),
                v         =as.double(v),
 		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		clusts    =as.integer(clusts),		
		iflag     =as.integer(iflag),
		imaxs     =as.integer(imaxs),
		imaxsc    =as.integer(imaxsc),
		imins     =as.integer(imins),
		iminsc    =as.integer(iminsc),
		index     =as.integer(index),				
		intcount  =as.integer(intcount),
		intcount2 =as.integer(intcount2),
		intind    =as.integer(intind),
		intind2   =as.integer(intind2),
		intposso  =as.integer(intposso),
		intpossn  =as.integer(intpossn),
		s         =as.integer(s),		
                seed      =as.integer(seed),
		betac     =as.double(betac),
		eta       =as.double(eta),
		etan      =as.double(etan),
		endp      =as.double(endp),
		endp2     =as.double(endp2),
		limr      =as.double(limr),		
                linfs     =as.double(linfs),
                linfsc    =as.double(linfsc),
                lsups     =as.double(lsups),
                lsupsc    =as.double(lsupsc),
		mass      =as.double(mass),
		prob      =as.double(prob),
                uvec      =as.double(uvec),
                vvec      =as.double(vvec),
		vnew      =as.double(vnew),
		wvec      =as.double(wvec),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		PACKAGE="DPpackage")	

         #########################################################################################
         # save state
         #########################################################################################
	
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+4))
 	 randsave<-matrix(foo$randsave,nrow=mcmc$nsave, ncol=(nrec+1))
 	 
 	 pnames<-NULL
 	 for(i in 1:p)
 	 {
 	    pnames<-c(pnames,namesxm[i+1])
 	 }
 	 
	 colnames(thetasave)<-c(pnames,"mu","sigma2","ncluster","alpha")

         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         
         colnames(randsave)<-qnames

	 model.name<-"Bayesian Semiparametric AFT Regression Model"
	 coeff<-rep(0,(p+4))
	 for(i in 1:(p+4)){
	 	coeff[i]<-mean(thetasave[,i])
	 }
	 names(coeff)<-c(pnames,"mu","sigma2","ncluster","alpha")
	
	 state <- list(alpha=foo$alpha,beta=foo$beta,v=foo$v,mu=foo$mu,sigma=foo$sigma)	
		      
	 save.state <- list(thetasave=thetasave,randsave=randsave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=foo$acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,cpo=foo$cpo,
	         nrec=nrec,p=p,possiP=possiP)
	
	 cat("\n\n")
	 class(z)<-c("DPsurvint")
	 z 
}


###                    
### Estimate the survival curve base on a fitted semiparametric aft model 
### for interval censored data.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 08-10-2006.


"predict.DPsurvint"<-
function(object,grid,xnew=NULL,hpd=TRUE, ...)
{
   if(is(object, "DPsurvint"))
   {
   	if(is.null(xnew))
   	{
   	     base<-1	
   	     npred<-1
   	     xnew<-rep(0,object$p)
   	     ngrid<-length(grid)
   	     nrec<-object$nrec
   	     nsave<-dim(object$save.state$thetasave)[1]
   	     mulog<-object$save.state$thetasave[,(object$p+1)]
	     sdlog<-sqrt(object$save.state$thetasave[,(object$p+2)])
	     alpha<-object$save.state$thetasave[,(object$p+4)]
	      
	     v<-matrix(object$save.state$randsave,nsave,nrec+1)
	     vpred<-v[,(nrec+1)]
	     v<-v[,1:nrec]
	     covn<-rep(0,1)
	     covn[1]<-"Baseline"
   	}
   	else
   	{
   	     base<-0
   	     npred<-dim(xnew)[1]
   	     pnew<-dim(xnew)[2]
   	     ngrid<-length(grid)
   	     nrec<-object$nrec
   	     nsave<-dim(object$save.state$thetasave)[1]
   	     
   	     mulog<-object$save.state$thetasave[,(object$p+1)]
   	     sdlog<-sqrt(object$save.state$thetasave[,(object$p+2)])
   	     alpha<- object$save.state$thetasave[,(object$p+4)]
   	     
   	     v <- matrix(object$save.state$randsave,nsave,nrec+1)
   	     vpred <- v[,(nrec+1)]
   	     v <- v[,1:nrec]

   	     if (object$p != pnew) 
   	     {
	       stop("Dimension of xnew is not the same that the design matrix
	        in the original model.\n")
	     }

	     covn <- rep(0,npred) 
	     
	     for(i in 1:npred)
	     {
			 covnw<-round(xnew[i,1],3)
                 
			 if(pnew>1)
			 {
				for(j in 2:pnew)
				{
					covnw<-paste(covnw,round(xnew[i,j],3),sep=";") 
				}
			}  
			covn[i]<-covnw  
	 	 }    
	 }    

 	 pm<-matrix(0,npred,ngrid)
	 pmed<-matrix(0,npred,ngrid)
	 psd<-matrix(0,npred,ngrid)
	 pstd<-matrix(0,npred,ngrid)
	 plinf<-matrix(0,npred,ngrid)
	 plsup<-matrix(0,npred,ngrid)
	     
	 expxb<-exp(xnew%*%t(object$save.state$thetasave[,1:object$p]))
	        
	 for(i in 1:npred)
	 {
	     surv <- matrix(0,nrow=ngrid,ncol=nsave)
	     for(j in 1:ngrid)
	     {
	         for(k in 1:nsave)
	         {
				 linf <- grid[j]*expxb[i,k]
				 mwork <- mulog[k]
				 swork <- sdlog[k]
				 awork <- alpha[k]
				 vwork <- v[k,]
				 Cdelta <- sum(vwork>linf)
				 Cparam <- awork*plnorm(linf,meanlog=mwork,sdlog=swork,
						  			    lower.tail=FALSE,log.p=FALSE)
				 
				 if(is.null(Cparam))
				 { 
					 cat(Cparam)
					 stop("Error")
				 }
				 if(is.null(Cdelta))
				 {
				     cat(Cdelta)
					 stop("Error") 
                 }    
				 
				 tmp1 <- (Cparam+Cdelta)
				 tmp2 <- alpha[k]+nrec-tmp1

				 tmp3 <- tmp1/(tmp1+tmp2)
				 if(tmp2<=0.1)
				 {
					 tmp3 <- 1
				 }
				 if(tmp1<=0.1)
				 {
					 tmp3 <- 0
				 }
				 if((tmp1 > 0.1) & (tmp2 > 0.1))
				 {
					 tmp3 <- rbeta(1,tmp1,tmp2) 
				 }
				 surv[j,k] <- tmp3
				 if(is.na(surv[j,k]))
				 {
				     cat(tmp1,"\n")
					 cat(tmp2,"\n")
					 stop("Error")
				 }	 
			 }
			 
			 pm[i,j]<-mean(surv[j,])
			 pmed[i,j]<-median(surv[j,])
			 psd[i,j]<-sqrt(var(surv[j,]))
			 pstd[i,j]<-sqrt(var(surv[j,]))/sqrt(nsave)
		  
			if(hpd==TRUE)
			{
				alow<-rep(0,2)
				aupp<-rep(0,2)
				sig<-0.05
				vec<-surv[j,]
				n<-length(vec)
				a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(sig),x=as.double(vec),
							alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
				plinf[i,j]<-a$alow[1]            
				plsup[i,j]<-a$aupp[1]            
			}
			else
			{
		        plinf[i,j]<-quantile(surv[j,],0.025)		
		        plsup[i,j]<-quantile(surv[j,],0.975)			                
			}
		}
   	 }
   	 
   	 dimnames(pm)<-list(covn,grid)
   	 dimnames(pmed)<-list(covn,grid)
   	 dimnames(psd)<-list(covn,grid)
   	 dimnames(pstd)<-list(covn,grid)
   	 dimnames(plinf)<-list(covn,grid)
   	 dimnames(plsup)<-list(covn,grid)
   	 
   	 out<-NULL
   	 out$pmean<-as.matrix(pm)
   	 out$pmedian<-as.matrix(pmed)
   	 out$psd<-as.matrix(psd)
   	 out$pstd<-as.matrix(pstd)
   	 out$plinf<-as.matrix(plinf)
   	 out$plsup<-as.matrix(plsup)
   	 out$xnew<-xnew
   	 out$vpred<-vpred
   	 out$grid<-grid
   	 out$npred<-npred
   	 out$base<-base
   	 out$covn<-covn
   	 class(out)<-c("predict.DPsurvint")
	 out
   }
}   


plot.predict.DPsurvint<-function(x,ask=TRUE,all=TRUE,band=FALSE,xlim=NULL,nfigr=1,nfigc=1, ...)
{

    if(is(x, "predict.DPsurvint"))
    {	
		if(x$base==1)
		{
		        par(ask = ask)
		        layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
			plot(density(x$vpred),xlim=xlim,xlab="Time",ylab="density",
			     main="Baseline Density")
			plot(c(0,x$grid),c(1,x$pmean[1,]),type="l",xlim=xlim,ylim=c(0,1),
			     xlab="Time",ylab="Survival",
			     main="Baseline Survival Curve",lty=1)
			lines(c(0,x$grid),c(1,x$plinf),lty=2)
			lines(c(0,x$grid),c(1,x$plsup),lty=2)
		}
		else
		{
			if(all==TRUE)
			{
			    plot(c(0,x$grid),c(1,x$pmean[1,]),type="l",xlim=xlim,ylim=c(0,1),
				 xlab="Time",ylab="Survival",
				 main="Survival Curves",lty=1)	
			    for(i in 2:x$npred)
			    {
			        lines(c(0,x$grid),c(1,x$pmean[i,]),lty=1)
			    }
			}
			else
			{
			    par(ask = ask)
			    layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
			    rnames<-rownames(x$pmean)
			    for(i in 1:x$npred)
		            {
		                title<-paste("Survival Curve for x=c(",rnames[i],sep="")	
		                title<-paste(title,")",sep="")	
		                plot(c(0,x$grid),c(1,x$pmean[i,]),type="l",xlim=xlim,ylim=c(0,1),
				     xlab="Time",ylab="Survival",
				     main=title,lty=1)
			        if(band==TRUE)lines(c(0,x$grid),c(1,x$plinf[i,]),lty=2)	  
			        if(band==TRUE)lines(c(0,x$grid),c(1,x$plsup[i,]),lty=2)	
		            }
			
			}
		}
    }		
}



###
### Tools for DPsurvint: anova, print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 15-12-2006.


"anova.DPsurvint"<-function(object, ...)
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
        paste("Response:", deparse(formula(object)[[2]]))), class = c("anovaPsCP",
        "data.frame"))
}



"print.DPsurvint"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\n\n")
    invisible(x)
}


"plot.DPsurvint"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
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


    if(is(x, "DPsurvint")){
        if(is.null(param))
	{
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:(n-1)){
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
	       {
	          hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
	       }
	       else
	       {
                  fancydensplot(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
               fancydensplot(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
            fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
        }
   }
}


"summary.DPsurvint"<-function(object, hpd=TRUE, ...) 
{
    dimen<-object$p
    coef.p<-object$coefficients[1:dimen]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])
    
    alpha<-0.05
    
    for(i in 1:dimen){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,i]))
        coef.m[i]<-median(object$save.state$thetasave[,i])
        vec<-object$save.state$thetasave[,i]
        n<-length(vec)
        
        if(hpd==TRUE)
        {
        
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                coef.l[i]<-a$alow[1]            
                coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    
    if(hpd==TRUE)
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

    ans$cpo<-object$cpo
    
    dimen<-2
    coef.p<-object$coefficients[(object$p+1):(object$p+dimen)]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    names(coef.l)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    names(coef.u)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    
    alpha<-0.05
    
    for(i in 1:dimen){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,(object$p+i)]))
        coef.m[i]<-median(object$save.state$thetasave[,(object$p+i)])
        vec<-object$save.state$thetasave[,(object$p+i)]
        n<-length(vec)

        if(hpd==TRUE)
        {
        
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                              alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                coef.l[i]<-a$alow[1]            
                coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd==TRUE)
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
    
    if(is.null(object$prior$a0))
    {
	dimen<-1	
    }
    else
    {
	dimen<-2	    
    }

    coef.p<-object$coefficients[(object$p+2+1):(object$p+2+dimen)]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[(object$p+2+1):(object$p+2+dimen)])
    names(coef.l)<-names(object$coefficients[(object$p+2+1):(object$p+2+dimen)])
    names(coef.u)<-names(object$coefficients[(object$p+2+1):(object$p+2+dimen)])
    for(i in 1:dimen){
         alow<-rep(0,2)
         aupp<-rep(0,2)
         coef.sd[i]<-sqrt(var(object$save.state$thetasave[,(object$p+2+i)]))
         coef.m[i]<-median(object$save.state$thetasave[,(object$p+2+i)])
         vec<-object$save.state$thetasave[,(object$p+2+i)]
         n<-length(vec)
        
         if(hpd==TRUE)
         {
             a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                              alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
             coef.l[i]<-a$alow[1]            
             coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd==TRUE)
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

    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryDPsurvint"
    return(ans)
}


"print.summaryDPsurvint"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
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


    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


