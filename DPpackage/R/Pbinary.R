### Pbinary.R
### Fit a parametric bernoulli regression model.
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 15-12-2006.
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

"Pbinary"<-
function(formula,link="logit",prior,mcmc,state,status,misc=NULL,
data=sys.frame(sys.parent()),na.action=na.fail) 
UseMethod("Pbinary")

"Pbinary.default"<-
function(formula,
         link="logit",
         prior,
         mcmc,
         state,
         status,
         misc=NULL,
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

  	 yobs<- model.response(mf,"numeric")
	 nrec<-length(yobs)
	 x<-as.matrix(model.matrix(formula))
	 p<-dim(x)[2]

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
         startp<-1+attr(tt, "intercept")
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
         # misclassification
         #########################################################################################
 
 	 if(is.null(misc))
	 {
	  	sens<-rep(1,nrec)
	  	spec<-rep(1,nrec)
	 }
	 else
	 {
 		sens<-misc$sens
	 	spec<-misc$spec
		if(length(sens)==1)sens<-rep(sens,nrec)
		if(length(spec)==1)spec<-rep(spec,nrec)
	 }

         #########################################################################################
         # mcmc specification
         #########################################################################################

         MLElogit<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-plogis(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }
         

         MLEprobit<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-pnorm(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }
 
         pcloglog<-function(x)
         {
            return(1-exp(-exp(x)))
         }


         MLEcloglog<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-pcloglog(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }


         MLEcauchy<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-pcauchy(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }

         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

          
         xmatrix<-x

         if(link=="logit")
         {
             linkn<-1
	     fit0<- MLElogit(xmatrix,yobs,sens,spec)
	 }
         if(link=="probit")
         {
             linkn<-2
	     fit0<- MLEprobit(xmatrix,yobs,sens,spec)
	 }
         if(link=="cloglog")
         {
             linkn<-3
	     fit0<- MLEprobit(xmatrix,yobs,sens,spec)
	 }
         if(link=="cauchy")
         {
             linkn<-4
	     fit0<- MLEcauchy(xmatrix,yobs,sens,spec)
	 }

         propv<- fit0$covb

         if(is.null(mcmc$tune))
         {
            tune=1
         }
         else
         {
            tune<-mcmc$tune
         }

         #########################################################################################
         # prior information
         #########################################################################################
	 
  	 betapm<-prior$beta0
	 betapv<-prior$Sbeta0
	
	 propv<-diag(mcmc$tune,p)%*%solve(solve(betapv)+solve(propv))%*%diag(mcmc$tune,p)

         #########################################################################################
         # parameters depending on status
         #########################################################################################

       	 if(status)
	 {
                beta<-fit0$beta
		eta <- x %*% beta
	 }	
      	 else
	 {
		beta<-state$beta
		eta <- x %*% beta
	 }
	 

         #########################################################################################
         # output
         #########################################################################################
         
	 thetasave <- matrix(0, nrow=mcmc$nsave, ncol=p)


         #########################################################################################
         # working space
         #########################################################################################

         acrate<-0	 
	 betac<-rep(0,p)
	 etan<-rep(0,nrec)
	 iflag<-rep(0,p)
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)
	 seed1<-sample(1:29000,1)
	 seed2<-sample(1:29000,1)
	 cpo<-rep(0,nrec)

         #########################################################################################
         # calling the fortran code
         #########################################################################################

 	 foo <- .Fortran("pbinary",
	 	link      =as.integer(linkn),
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		sens      =as.double(sens),
		spec      =as.double(spec),
		x         =as.double(x),
		yobs      =as.integer(yobs),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		acrate    =as.double(acrate),
		thetasave =as.double(thetasave),		
		cpo       =as.double(cpo),		
		beta      =as.double(beta),
		betac     =as.double(betac),
		eta       =as.double(eta),
		etan      =as.double(etan),
		iflag     =as.integer(iflag),
		seed1     =as.integer(seed1),
		seed2     =as.integer(seed2),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		PACKAGE="DPpackage")	


         #########################################################################################
         # save state
         #########################################################################################
	
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p))

         model.name<-"Bayesian parametric binary regression model"		

  	 colnames(thetasave)<-c(dimnames(x)[[2]])
  	 
  	 coeff<-rep(0,p)
	 for(i in 1:p){
	     coeff[i]<-mean(thetasave[,i])
	 }  
	 names(coeff)<-c(dimnames(x)[[2]])
	
	
	 state <- list(beta=foo$beta)
				  
	 save.state <- list(thetasave=thetasave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=foo$acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
	         cpo=foo$cpo,p=p,link=link,x=x,possiP=possiP)
	         
	 cat("\n\n")

	 class(z)<-c("Pbinary")
	 return(z) 
}


###                    
### Estimate the probability curve for a fitted parametric binary 
### regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 01-07-2006.

predict.Pbinary<-function(object,xnew=NULL,hpd=TRUE, ...)
{
   pcloglog<-function(x)
   {
      return(1-exp(-exp(x)))
   }
   
   sd<-function(x)
   {
     return(sqrt(var(x)))
   }

   std<-function(x)
   {
     n<-length(x)
     return(sqrt(var(x))/sqrt(n))
   }

   hpdf<-function(x)
   {
      alow<-rep(0,2)
      aupp<-rep(0,2)
      sig<-0.05
      n<-length(x)
      a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(sig),x=as.double(x),
                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
      return(c(a$alow[1],a$aupp[1]))
   }
   
   nhpdf<-function(x)
   {
      return(c(quantile(x,0.025),quantile(x,0.975)))
   }
   
   
   if(is.null(xnew))
   {
      xnew<-object$x
   }

   if(is(object, "Pbinary"))
   {
   	 npred<-dim(xnew)[1]
   	 pnew<-dim(xnew)[2]
   	 nrec<-object$nrec
   	 if (object$p != pnew) 
   	 {
	   stop("Dimension of xnew is not the same that the design matrix
	    in the original model.\n")
	 }

	 plinf<-rep(0,npred)
	 plsup<-rep(0,npred)
	     
	 lp<-xnew%*%t(object$save.state$thetasave[,1:object$p])

         if(object$link=="logit")
         {
	     prob<-plogis(lp) 
	 }
         if(object$link=="probit")
         {
	     prob<-pnorm(lp) 
	 }
         if(object$link=="cloglog")
         {
	     prob<-pcloglog(lp) 
	 }
         if(object$link=="cauchy")
         {
	     prob<-pcauchy(lp) 
	 }
         
         pm<-apply(prob,1,mean)
         pmed<-apply(prob,1,median)
         psd<-apply(prob,1,sd)
         pstd<-apply(prob,1,std)
         
         if(hpd=="TRUE")
         {
            phpd<-apply(prob,1,hpdf)
         }
         else
         {
            phpd<-apply(prob,1,nhpdf)         
         }

	 covn<-rep(0,npred) 
	 
	 for(i in 1:npred)
	 {
             covnw<-round(xnew[i,1],3)
             for(j in 2:pnew){
                 covnw<-paste(covnw,round(xnew[i,j],3),sep=";") 
             }
             covn[i]<-covnw  
         }    

         plinf<-phpd[1,] 
         plsup<-phpd[2,]

   	 names(pm)<-covn
   	 names(pmed)<-covn
   	 names(psd)<-covn
   	 names(pstd)<-covn
   	 names(plinf)<-covn
   	 names(plsup)<-covn
   	 
   	 out<-NULL
   	 out$pmean<-pm
   	 out$pmedian<-pmed
   	 out$psd<-psd
   	 out$pstd<-pstd
    	 out$plinf<-plinf
    	 out$plsup<-plsup
         return(out)
   }	    
}


###
### Tools for Pbinary: anova, print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 15-12-2006.


"anova.Pbinary"<-function(object, ...)
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

	

"print.Pbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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



"summary.Pbinary"<-function(object, hpd=TRUE, ...) 
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
    
    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryPbinary"
    return(ans)
}



"print.summaryPbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}




"plot.Pbinary"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
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


    if(is(x, "Pbinary")){
        if(is.null(param))
	{
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:n){
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






