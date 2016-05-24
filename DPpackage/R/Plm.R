### Plm.R                   
### Fit a parametric regression model.
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 15-03-2009.
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

"Plm"<-
function(formula,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("Plm")

"Plm.default"<-
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
         # data and model structure
         #########################################################################################

 	 y <-  model.response(mf,"numeric")
  	 nrec <- length(y)
  	 x <- model.matrix(formula)
  	 p <- ncol(x)

  	 nfixed <- p
 
         #########################################################################################
         # Elements for Pseudo Countour Probabilities' computation
         #########################################################################################
         Terms <- if (missing(data)) 
              terms(formula)
         else terms(formula, data = data)

         possiP <- NULL
         if(nfixed>0)
         {
            mat <- attr(Terms,"factors")
            namfact <- colnames(mat)
            nvar <- dim(mat)[1]
            nfact <- dim(mat)[2]
            possiP <- matrix(0,ncol=2,nrow=nfact)
            if (missing(data)) dataF <- model.frame(formula=formula,xlev=NULL)
               dataF <- model.frame(formula=formula,data,xlev=NULL)
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

            startp<-1+1

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

         betapm <- prior$beta0
         betapv <- prior$Sbeta0

         sbeta0i <- solve(betapv)
         ms <- sbeta0i%*%betapm

         tau1<-prior$tau1
         tau2<-prior$tau2

         #########################################################################################
         # mcmc specification
         #########################################################################################

         if(missing(mcmc))
         {
            nburn <- 1000
            nsave <- 1000
            nskip <- 0
            ndisplay <- 100
            mcmcvec<-c(nburn,nskip,ndisplay)
         }
         else
         {
            mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
            nsave<-mcmc$nsave
         }

         #########################################################################################
         # output
         #########################################################################################
         
         cpo <- matrix(0,nrow=nrec,ncol=2)
	 thetasave <- matrix(0, nrow=nsave, ncol=p+2)

         #########################################################################################
         # parameters depending on status
         #########################################################################################

	 if(status==TRUE)
	 {
            fit0 <- lm(formula)
            beta <- coefficients(fit0)
            sigma2 <- sum((resid(fit0))**2)/(nrec-p)
 	 }	
	 if(status==FALSE)
	 {
	    beta <- state$beta
	    sigma2 <- state$sigma2
	 }

         #########################################################################################
         # working space
         #########################################################################################
         betam <- rep(0,p)
         seed <- c(sample(1:29000,1),sample(1:29000,1))
 	 iflagp <- rep(0,p)
	 workmh <- rep(0,(p*(p+1)/2))
	 xtx <- matrix(0,nrow=p,ncol=p)
	 xty <- rep(0,p)
	 
         #########################################################################################
         # calling the fortran code
         #########################################################################################

           foo <- .Fortran("plm",
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		x         =as.double(x),	 	
		y         =as.double(y),
                tau1      =as.double(tau1),
                tau2      =as.double(tau2),
                ms        =as.double(ms),
                sbeta0i   =as.double(sbeta0i),
                beta      =as.double(beta), 
                sigma2    =as.double(sigma2),
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
                cpo       =as.double(cpo),
		thetasave =as.double(thetasave),
		iflagp    =as.integer(iflagp),
		betapm    =as.double(betapm),		
		workmh    =as.double(workmh),
		xtx       =as.double(xtx),
		xty       =as.double(xty),
                seed      =as.integer(seed),
		PACKAGE="DPpackage")	                  
         
         #########################################################################################
         # save state
         #########################################################################################

         cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo<-cpom[,1]         
         fso<-cpom[,2]
	
 	 thetasave <- matrix(foo$thetasave,nrow=nsave, ncol=(p+2))
 	 
	 colnames(thetasave)<-c(dimnames(x)[[2]],"sigma2","LPML")

	 model.name<-"Bayesian Linear Regression Model"

         coeff<-apply(thetasave,2,mean)[1:(p+1)]		
	
	 state <- list(beta=foo$beta,
                       sigma2=foo$sigma2)	
		      
	 save.state <- list(thetasave=thetasave)


	 z<-list(modelname=model.name,
                 coefficients=coeff,
                 call=cl,
	         prior=prior,
                 mcmc=mcmc,
                 state=state,
                 save.state=save.state,
                 cpo=cpo,
                 fso=fso,
	         nrec=nrec,
                 p=p,
                 possiP=possiP)
	
	 cat("\n\n")
	 class(z)<-c("Plm")
	 z 
}


###
### Tools for Plm: anova, print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2007
### Last modification: 14-08-2007.


"anova.Plm"<-function(object, ...)
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




"print.Plm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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
    cat("\nNumber of Preditors:",x$p)
    cat("\n\n")
    invisible(x)
}


"plot.Plm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
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


    if(is(x, "Plm")){
        if(is.null(param))
	{
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:n)
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
            par(ask = ask)
            layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
            title1<-paste("Trace of",pnames[poss],sep=" ")
            title2<-paste("Density of",pnames[poss],sep=" ")       
            plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
            fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
        }
   }
}


"summary.Plm"<-function(object, hpd=TRUE, ...) 
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

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo <- object$cpo

### Regression

    dimen1 <- object$p+1
    mat<- thetasave[,1:dimen1] 

    coef.p <- object$coefficients[1:dimen1]
    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

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
    ans$coefficients <- coef.table
    ans$p <- object$p
    ans$cpo <- object$cpo
    ans$nrec <- object$nrec

    class(ans) <- "summaryPlm"
    return(ans)
}



"print.summaryPlm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) 
    {
        cat("\nRegression coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$prec)) 
    {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat(" \n")

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Preditors:",x$p)
    cat("\n\n")
    invisible(x)
}


