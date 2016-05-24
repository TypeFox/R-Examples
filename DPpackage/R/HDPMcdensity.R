### HDPMcdensity.R                   
### Fit a hierarchical mixture of DPM of normals model for conditional density 
### estimation
###
### Copyright: Alejandro Jara and Peter Mueller, 2008-2012.
###
### Last modification: 07-06-2008.
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
###      Peter Mueller
###      Department of Mathematics
###      The University of Texas Austin
###      1, University Station, C1200 
###      Austin TX 78712, USA
###      Voice: (512) 471-7168  URL  : http://www.math.utexas.edu/users/pmueller
###      Fax  : (512) 471-9038  Email: pmueller@math.utexas.edu
###

"HDPMcdensity"<-function(formula,study,xpred,ngrid=100,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("HDPMcdensity")

"HDPMcdensity.default"<-
function(formula,
         study,
         xpred,
  		 ngrid=100,
         prior,
         mcmc,
         state,
         status,
         data=sys.frame(sys.parent()),
         na.action=na.fail,
         work.dir=NULL)  
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

         y <- as.matrix(model.response(model.frame(formula, data=data)))
         nrec <- nrow(y)
		 nvar <- ncol(y)
         x <- as.matrix(model.matrix(formula, data=data))
         p <- ncol(x)
         p <- p-1

		 m <- mcall <- cl <- match.call()
         nm <- names(m)[-1]
         keep <- is.element(nm, c("data", "na.action"))
         for (i in nm[!keep]) m[[i]] <- NULL

         allvars <- c(all.vars(formula), all.vars(study))

         Terms <- if (missing(data)) 
               terms(formula)
		 else terms(formula, data = data)

         cl$fixed <- eval(formula)
         cl$study <- eval(study)
         m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
         environment(m$formula) <- environment(formula)
         m$drop.unused.levels <- TRUE
         m[[1]] <- as.name("model.frame")
         mf <- eval.parent(m)

	   #########################################################################################
       # data structure
       #########################################################################################
         mf <- as.matrix(mf)
         nrec <- nrow(mf)
         y <- mf[,1:nvar]
         x <- mf[,(nvar+1):(nvar+p)]
		 study <- mf[,(nvar+p+1)]
		 namesxm <- colnames(x)
         ntvar <- nvar + p

         nstudy <- max(study)

         if (sum(abs(study-round(study)))>0)
         { # not integers
             cat("\n *** Error: study indicators need to be integers 1...nstudies.\n")
             return(-1)
         }
         if (length(unique(study)) != max(study))
         {
            cat("\n *** Error: studies need to be indexed 1...nstudies.\n")
            return(-1)
         }

       #########################################################################################
       # prediction
       #########################################################################################
         npred <- nrow(xpred)
         if ( (ncol(xpred) != p) )
		 {
              cat("\n *** Error: dim(xpred) != (npred x p).\n")
              return(-1)
         }

       #########################################################################################
       # prior specification
       #########################################################################################
         if(is.null(prior$pe1))
         { 
            pe1 <- 0.1
         }
         else
         {
            pe1 <- prior$pe1
         } 

         if(is.null(prior$pe0))
         { 
            pe0 <- 0.1
         }
         else
         {
            pe0 <- prior$pe0
         } 

         pieps <- c(pe0,pe1)

         if (pe1+pe0 >= 1.0)
         {
             cat("\n *** Error: need pe0+pe1 < 1.\n")
             return(-1)
         }

         if(is.null(prior$eps))
         { 
            sameps <- 1
            ae <- prior$ae
            be <- prior$be
            if ( (ae <= 0) | (be <= 0))
            {
             cat("\n *** Error: need ae > 0 and be >0.\n")
             return(-1)
            }
            aebe <- c(ae,be)
            eps <- runif(1)
         }
         else
         {
            sameps <- 0
            aebe <- rep(0,2)
            eps <- prior$eps
         } 

		 if(is.null(prior$alpha))
         {
            if(is.null(prior$a))
            { 
               a0 <- rep(1,nstudy+1)
            }
            else
            {
               a0 <- prior$a0
               if(length(a0)!=(nstudy+1))
               {
                  cat("\n *** Error: length of a0 must be nstudies+1.\n")
                  return(-1)
               }
            } 
         
            if(is.null(prior$b0))
            { 
               b0 <- rep(1,nstudy+1)
            }
            else
            {
               b0 <- prior$b0
               if(length(b0)!=(nstudy+1))
               {
                  cat("\n *** Error: length of b0 must be nstudies+1.\n")
                  return(-1)
               }
            }
            alpha <- rgamma(nstudy+1,shape=a0,rate=b0)
            a0b0 <- cbind(a0,b0)
         }
         else
         {
            alpha <- prior$alpha
            if(length(alpha)!=(nstudy+1))
            {
			   cat("\n *** Error: length of alpha must be nstudies+1.\n")
               return(-1)
            }
            a0b0 <- matrix(-1,nrow=nstudy+1,ncol=2)
         }


		 if(is.null(prior$sigma))
         {
			nu <- prior$nu
            if (nu < ntvar+2)
            {    
                cat(" *** Warning: should use nu > nvar+p+2 for good mixing MCMC.\n")
            }
            tinv <- prior$tinv
            if ( (nrow(tinv) != ntvar) | (ncol(tinv) != ntvar))
            {
               cat("\n *** Error: dim(tinv) != ((nvar + p) x (nvar + p)).\n")
               cat(" Note: dim(tinv) includes the covariates!.\n")
               return(-1)
            }
			sigma <- var(cbind(y,x))
         }
         else
         {
            sigma <- prior$sigma
            if ( (nrow(sigma) != ntvar) | (ncol(sigma) != ntvar))
            {
               cat("\n *** Error: dim(sigma) != ((nvar+p) x (nvar+p)).\n")
               cat(" Note: dim(sigma) includes the covariates!.\n")
               return(-1)
            }
            nu <- -1
            tinv <- diag(1,ntvar)
		 }

         if(is.null(prior$mub))
         {
            sammu <- 1
            mub <- apply(cbind(y,x),2,mean)
            m0 <- prior$m0
            prec0 <- solve(prior$S0)

            if(length(m0) != ntvar)
            {
               cat("\n *** Error: length(m0) != (nvar+p).\n")
               cat(" Note: length(m0) includes the covariates!.\n")
               return(-1)
            }

            if ( (nrow(prec0) != ntvar) | (ncol(prec0) != ntvar))
            {
               cat("\n *** Error: dim(S0) != ((nvar+p) x (nvar+p)).\n")
               cat(" Note: dim(S0) includes the covariates!.\n")
               return(-1)
            }
         }
         else
         {
            sammu <- 0
            mub <- prior$mub
            if(length(mub) != ntvar)
            {
               cat("\n *** Error: length(mub) != (nvar+p).\n")
               cat(" Note: length(mub) includes the covariates!.\n")
               return(-1)
            }
			m0 <- rep(0,ntvar)
		    prec0 <- diag(1,ntvar)
         }


         if(is.null(prior$sigmab))
         {
			samsb <- 1
            nub <- prior$nub
            if (nub < ntvar+2)
            {    
                cat(" *** Warning: should use nub > nvar+p+2 for good mixing MCMC.\n")
            }
            tbinv <- prior$tbinv
            if ( (nrow(tbinv) != ntvar) | (ncol(tbinv) != ntvar))
            {
               cat("\n *** Error: dim(tbinv) != ((nvar+p) x (nvar+p)).\n")
               cat(" Note: dim(tbinv) includes the covariates!.\n")
               return(-1)
            }
			sigmab <- var(cbind(y,x))
		 }
         else
         { 
            samsb <- 0
            sigmab <- prior$sigmab
            nub <- 0 
            tbinv <- diag(1,ntvar)
            if ( (nrow(sigmab) != ntvar) | (ncol(sigmab) != ntvar))
            {
               cat("\n *** Error: dim(sigmab) != ((nvar+p) x (nvar+p)).\n")
               cat(" Note: dim(sigmab) includes the covariates!.\n")
               return(-1)
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
         densmc <- matrix(0,nrow=(nstudy+1),ncol=npred*nvar*ngrid)
         densms <- matrix(0,nrow=nstudy,ncol=npred*nvar*ngrid)
         ntotal <- ntvar*(ntvar+1)/2+1+nstudy+1+ntvar+ntvar*(ntvar+1)/2
         thetasave <- matrix(0,nrow=nsave,ncol=ntotal) 
         randsave <- matrix(0,nrow=nsave,ncol=nstudy)

       #########################################################################################
       # parameters depending on status
       #########################################################################################
        
    	 if(status==TRUE)
	     { 
            ncluster <- nstudy
            ss <- study
            sc <- rep(0,nrec)
            sc[1:nstudy] <- seq(1,nstudy)
            muclus <- matrix(0,nrow=(nrec+100),ncol=ntvar)
            for(i in 1:ncluster)
            {
                muclus[i,1:ntvar] <- apply(cbind(y[study==i,],x[study==i,]),2,mean)
            }
         }
	 
      	 if(status==FALSE)
	     {
            ncluster <- state$ncluster
            ss <- state$ss
            sc <- state$sc
            alpha <- state$alpha
            muclus <- state$muclus
            sigma <- state$sigma
            mub <- state$mub
            sigmab <- state$sigmab
            eps <- state$eps
	    }    

       #########################################################################################
       # working space
       #########################################################################################
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)

         cstrt <- matrix(0,nrow=nrec,ncol=nrec)
         ccluster <- rep(0,nrec)
         scstrt <- matrix(0,nrow=nstudy+1,ncol=nrec)
         sccluster <- rep(0,nstudy+1)
         iflagp <- rep(0,ntvar)
         ywork <- rep(0,ntvar)
         muwork <- rep(0,ntvar)
         sigmainv <- matrix(0,nrow=ntvar,ncol=ntvar)
		 sigmabinv <- matrix(0,nrow=ntvar,ncol=ntvar)
         sigmawork <- matrix(0,nrow=ntvar,ncol=ntvar)
         prob <- rep(0,nrec+100)
         quadf <- matrix(0,nrow=ntvar,ncol=ntvar)
         workm1 <- matrix(0,nrow=ntvar,ncol=ntvar)
         workmh1 <- rep(0,ntvar*(ntvar+1)/2)
         workmh2 <- rep(0,ntvar*(ntvar+1)/2)
         workv1 <- rep(0,ntvar)

         iflagx <- rep(0,p)
         densw <- rep(0,nvar*ngrid)
         mubar <- rep(0,nvar)
         sigmabar <- matrix(0,nrow=nvar,ncol=nvar)
         sigmabarb <- matrix(0,nrow=nvar,ncol=nvar)
         workmx1 <- matrix(0,nrow=p,ncol=p)
         workmx2 <- matrix(0,nrow=nvar,ncol=p)
         workmx3 <- matrix(0,nrow=nvar,ncol=p)
         workmx3b <- matrix(0,nrow=nvar,ncol=p)
         workvx1 <- rep(0,p)
         workvx2 <- rep(0,nvar)
         workvx3 <- rep(0,p)

       #########################################################################################
       # grid for prediction
       #########################################################################################

         grid <- NULL
         for(i in 1:nvar)
         {
             left <- min(y[,i]) - 0.5*sqrt(var(y[,i]))
			 right <- max(y[,i]) + 0.5*sqrt(var(y[,i]))
             tmp <- seq(left,right,len=ngrid)
             grid <- cbind(grid,tmp)
		 }

       #########################################################################################
       # change working directory (if requested..)
       #########################################################################################
         if(!is.null(work.dir))
         {
            cat("\n Changing working directory to ",work.dir,"\n")
            old.dir <- getwd()  # by default work in current working directory
            setwd(work.dir)
         }

       #########################################################################################
       # calling the fortran code
       #########################################################################################
         mcmcvec <- c(mcmcvec,sameps,sammu,samsb)

		 foo <- .Fortran(  "hdpmcdensity",
							nrec       =as.integer(nrec),
							nvar       =as.integer(nvar),
							p		   =as.integer(p),
							ntvar      =as.integer(ntvar),
							nstudy	   =as.integer(nstudy),	
							study      =as.integer(study),
							y		   =as.double(y),
							x		   =as.double(x),
							ngrid      =as.integer(ngrid),
							npred      =as.integer(npred),
							grid       =as.double(grid),
							xpred      =as.double(xpred),
							nu		   =as.integer(nu),
							nub		   =as.integer(nub),
							a0b0	   =as.double(a0b0),
							tinv	   =as.double(tinv),
							m0		   =as.double(m0),
							prec0	   =as.double(prec0),	
							tbinv	   =as.double(tbinv),
							aebe	   =as.double(aebe),	
							pieps	   =as.double(pieps),
							ncluster   =as.integer(ncluster),
							ss		   =as.integer(ss),
							sc		   =as.integer(sc),
							alpha	   =as.double(alpha),
							muclus     =as.double(muclus),
							sigma      =as.double(sigma),
							mub		   =as.double(mub),
							sigmab	   =as.double(sigmab),
							eps		   =as.double(eps),
							mcmc       =as.integer(mcmcvec),
							nsave	   =as.integer(nsave),
                            seed       =as.integer(seed),
							densmc	   =as.double(densmc),	
							densms	   =as.double(densms),
							thetasave  =as.double(thetasave),
							randsave   =as.double(randsave),
							cstrt	   =as.integer(cstrt),
							ccluster   =as.integer(ccluster),	
							scstrt	   =as.integer(scstrt),	 	
							sccluster  =as.integer(sccluster),
							iflagp	   =as.integer(iflagp),
							ywork      =as.double(ywork),
							muwork     =as.double(muwork),
							sigmainv   =as.double(sigmainv),
							sigmabinv  =as.double(sigmabinv),
							sigmawork  =as.double(sigmawork),	
							prob	   =as.double(prob),
							quadf	   =as.double(quadf),
							workm1	   =as.double(workm1),
							workmh1	   =as.double(workmh1),
							workmh2	   =as.double(workmh2),
							workv1	   =as.double(workv1),
							iflagx	   =as.integer(iflagx),
							densw	   =as.double(densw),
							mubar	   =as.double(mubar),
							sigmabar   =as.double(sigmabar),
							sigmabarb  =as.double(sigmabarb),
							workmx1	   =as.double(workmx1),
							workmx2	   =as.double(workmx2),
							workmx3	   =as.double(workmx3),
							workmx3b   =as.double(workmx3b),
							workvx1    =as.double(workvx1),
							workvx2	   =as.double(workvx2),
                            workvx3	   =as.double(workvx3),
		                    PACKAGE    ="DPpackage")

       #########################################################################################
       # save state
       #########################################################################################

         varnames <- colnames(cbind(y,x))

         model.name <- "Hierarchical Mixture of DPM of normals model"		

         densmc <- matrix(foo$densmc,nrow=(nstudy+1),ncol=npred*nvar*ngrid)
         densms <- matrix(foo$densms,nrow=nstudy,ncol=npred*nvar*ngrid)
         ntotal <- ntvar*(ntvar+1)/2+1+nstudy+1+ntvar+ntvar*(ntvar+1)/2
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=ntotal) 
         randsave <- matrix(foo$randsave,nrow=nsave,ncol=nstudy)
         
         state <- list(	ncluster=foo$ncluster,
						ss=foo$ss,
						sc=foo$sc,
						alpha=foo$alpha,
						muclus=matrix(foo$muclus,nrow=(nrec+100),ncol=ntvar),
						sigma=matrix(foo$sigma,nrow=ntvar,ncol=ntvar),
						mub=foo$mub,
						sigmab=matrix(foo$sigmab,nrow=ntvar,ncol=ntvar),
						eps=foo$eps)

         pnames <- NULL
		 for(i in 1:ntvar)
		 {
			 for(j in i:ntvar)
			 {
				 pnames <- c(pnames,paste("sigma",i,j,sep=""))
			 }
		 }
         pnames <- c(pnames,"eps")

         for(i in 1:nstudy)
         {
             pnames <- c(pnames,paste("alpha",i,sep=""))
         }
         pnames <- c(pnames,"alpha0")


         for(i in 1:ntvar)
         {
             pnames <- c(pnames,paste("mub",i,sep=""))
         }
		 for(i in 1:ntvar)
		 {
			 for(j in i:ntvar)
			 {
				 pnames <- c(pnames,paste("sigmab",i,j,sep=""))
			 }
		 }
         colnames(thetasave) <- pnames

		 save.state <- list(thetasave=thetasave,
							randsave=randsave)
  
         coeff <- apply(thetasave,2,mean)


       # return error code
         if(!is.null(work.dir))
         {
           cat("\n Changing working directory back to ",old.dir,"\n")
           setwd(old.dir)
         }


	     z <- list(modelname=model.name,
	               call=cl,
				   prior=prior,
				   mcmc=mcmc,
				   state=state,
                   save.state=save.state,
                   coefficients=coeff,
                   work.dir=work.dir,
                   nrec=nrec,
                   nvar=nvar,
                   p=p,
                   npred=npred,
				   xpred=xpred,
				   ntvar=ntvar,
				   densmc=densmc,
				   densms=densms,
				   grid=grid,
				   ngrid=ngrid,
                   nstudy=nstudy,
				   varnames=varnames)
                 
         cat("\n\n")        

         class(z)<-c("HDPMcdensity")
         return(z) 
}


###                    
### Tools
###
### Copyright: Alejandro Jara, 2009
### Last modification: 25-09-2009.
###


"print.HDPMcdensity" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nNumber of subjects:",x$nrec)
    cat("\nNumber of variables:",x$nvar)
    cat("\nNumber of predictors:",x$p)
    cat("\nNumber of studies:",x$nstudy)
    cat("\n\n")
    invisible(x)
}



"summary.HDPMcdensity"<-function(object, hpd=TRUE, ...) 
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

### kernel parameters

    ntvar <- object$ntvar 
    dimen1 <- ntvar*(ntvar+1)/2

    if(dimen1>1)
    {
       mat <- thetasave[,1:dimen1]
    }
    else
    {
	   mat <- matrix(thetasave[,1:dimen1],ncol=1)
    }

    coef.p <- object$coefficients[1:dimen1]
    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

    if(hpd){             
         limm <- apply(mat, 2, hpdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }
    else
    {
         limm <- apply(mat, 2, pdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
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

    ans$kernel <- coef.table

### eps

	dimen2 <- 1
	mat <- matrix(thetasave[,dimen1+1],ncol=1)

    coef.p <- object$coefficients[dimen1+1]
    coef.m <- apply(mat, 2, median)    
    coef.sd <- apply(mat, 2, sd)
    coef.se <- apply(mat, 2, stde)

    if(hpd){             
         limm <- apply(mat, 2, hpdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }
    else
    {
         limm <- apply(mat, 2, pdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
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

    ans$eps <- coef.table


### Precision parameter

    nstudy <- object$nstudy 
	dimen3 <- nstudy+1
    coef.p <- object$coefficients[(dimen1+dimen2+1):(dimen1+dimen2+dimen3)] 
	mat <- thetasave[,(dimen1+dimen2+1):(dimen1+dimen2+dimen3)]

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

	ans$prec <- coef.table

### Baseline Information

    dimen4 <- ntvar+ntvar*(ntvar+1)/2

    mat <- thetasave[,(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)]

	coef.p <- object$coefficients[(dimen1+dimen2+dimen3+1):(dimen1+dimen2+dimen3+dimen4)]
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


    ans$nrec <- object$nrec
    ans$nvar <- object$nvar
    ans$p <- object$p
    ans$nstudy <- object$nstudy

    class(ans) <- "summaryHDPMcdensity"
    return(ans)
}


"print.summaryHDPMcdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    if (length(x$kernel)) {
        cat("\nNormal kernel covariance:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

	if (length(x$eps)) {
        cat("\nWeight parameter:\n")
        print.default(format(x$eps, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

	if (length(x$prec)) {
        cat("\nPrecision parameters:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$base)) {
        cat("\nBaseline distribution:\n")
        print.default(format(x$base, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nNumber of subjects:",x$nrec)
    cat("\nNumber of variables:",x$nvar)
    cat("\nNumber of predictors:",x$p)
    cat("\nNumber of studies:",x$nstudy,"\n")
    cat("\n\n")
    invisible(x)
}



"plot.HDPMcdensity" <- function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "HDPMcdensity"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:n)
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
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

        }
		else
        {
			coef.p <- x$coefficients
			n <- length(coef.p)
			pnames <- names(coef.p)
			poss <- 0 
			for(i in 1:n)
			{
				if(pnames[i]==param)poss <- i
			}
            if(poss==0 && param !="prediction")             
			{
				stop("This parameter is not present in the original model.\n")
			}
	    
			par(ask = ask)
			layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))

			if(param !="prediction")
			{
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
            else
		    {


            }                

        }
   }

}


predict.HDPMcdensity <- function(object,pred,i,r,ask=TRUE,nfigr=2,nfigc=2, ...)
{
	   nvar <- object$nvar
       ntvar <- object$ntvar
       p <- object$p
       npred <- object$npred
       xpred <- object$xpred
       ngrid <- object$ngrid
       nstudy <- object$nstudy
	   grid <- object$grid
       nn <- object$varnames

       par(ask = ask)
       layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))

       if(i==0)
       {
          start <- ngrid*nvar*(pred-1)+1
          end <- ngrid*nvar*(pred-1)+ngrid
          for(j in 1:nvar)
          {
              tmp1 <- paste("f","0",sep="")
              tmp2 <- paste("Variable",nn[j],sep="=")
              tit <- paste(tmp1,tmp2,sep=" ; ")
              plot(grid[,j],object$densmc[nstudy+1,start:end],type="l",main=tit,xlab="values",ylab="density")
              start <- end+1
              end   <- end+ngrid 
         }
       }
       else
	   {
    
	     if(r==1)
         {
            start <- ngrid*nvar*(pred-1)+1
            end <- ngrid*nvar*(pred-1)+ngrid
            for(j in 1:nvar)
            {
                tmp1 <- paste("f",i,sep="")
                tmp2 <- paste("Variable",nn[j],sep="=")
                tit <- paste(tmp1,tmp2,sep=" ; ")
                plot(grid[,j],object$densms[i,start:end],type="l",main=tit,xlab="values",ylab="density")
                start <- end+1
                end   <- end+ngrid 
           }
        }
        else
        {
            start <- ngrid*nvar*(pred-1)+1
            end <- ngrid*nvar*(pred-1)+ngrid
            for(j in 1:nvar)
            {
				tmp1 <- paste("f",i,sep="")
                tmp2 <- paste("Variable",nn[j],sep="=")
                tit <- paste(tmp1,tmp2,sep=" ; ")
                plot(grid[,j],object$densmc[i,start:end],type="l",main=tit,xlab="values",ylab="density")
                start <- end+1
                end   <- end+ngrid 
           }
        }
     }
}
