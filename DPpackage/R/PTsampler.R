### PTsampler.R                   
### Generates a MCMC sampler from an user-specified unormilized multivariate
### density.
###
### Copyright: Alejandro Jara, 2009-2012.
###
### Last modification: 16-12-2011.
###
### In this version an error detected by Davor Cubranic has been corrected.
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


PTsampler <- function(ltarget,dim.theta,mcmc=NULL,support=NULL,pts.options=NULL,status=TRUE,state=NULL)

UseMethod("PTsampler")

PTsampler.default <- function(ltarget,dim.theta,mcmc=NULL,support=NULL,pts.options=NULL,status=TRUE,state=NULL)
{

  ######################################################################################
  # Internal functions
  ######################################################################################

	rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
						 method = c("eigen", "svd", "chol")) 
	{
		if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
			stop("sigma must be a symmetric matrix")
		}
		if (length(mean) != nrow(sigma)) {
			stop("mean and sigma have non-conforming size")
		}
		sigma1 <- sigma
		dimnames(sigma1) <- NULL
		if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
			warning("sigma is numerically not symmetric")
		}
		method <- match.arg(method)
		if (method == "eigen") {
			ev <- eigen(sigma, symmetric = TRUE)
			if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
				warning("sigma is numerically not positive definite")
			}
			retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
			t(ev$vectors)
		}
		else if (method == "svd") {
			sigsvd <- svd(sigma)
			if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
				warning("sigma is numerically not positive definite")
			}
			retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
		}
		else if (method == "chol") {
			retval <- chol(sigma, pivot = TRUE)
			o <- order(attr(retval, "pivot"))
			retval <- retval[, o]
		}
		retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
		retval <- sweep(retval, 2, mean, "+")
		retval
	}

	PTS.warmup <- function(aa) function(t)
	{
		#################################
		# Step 0	
		#################################
		  x <- aa$x
		  L1 <- aa$L1
		  cpar <- aa$cpar
		  u <- aa$u
   		  uinv <- aa$uinv
		#################################
		# Step 1	
		#################################
		  npoints <- nrow(x)
		  nvar <- ncol(x)
		
		  m <- apply(x,2,mean)
		  S <- var(x)

		  foo1 <- eigen(S)
		  sqrtLambda <- diag(sqrt(foo1$values))
		  nsqrtLambda <- diag(1/sqrt(foo1$values))
		  M <- foo1$vectors
	
		#################################
		# Step 2	
		#################################
		  A <- matrix(rnorm(nvar*nvar),nrow=nvar,ncol=nvar)
		  out <- qr(A)
		  Q <- qr.Q(out)
		  R <- qr.R(out)
          O <- Q%*%diag(sign(diag(R)))

		  u <- M%*%sqrtLambda%*%O
	 	  uinv <- t(M%*%nsqrtLambda%*%O)
			  
		  std.points <- function(x)
		  { 
			  uinv%*%(x-m)
		  }	  
		  z <- t(apply(x,1,std.points))
			  
		#################################
		# Step 3	
		#################################
		  ntint <- 2**(nvar*nlevel)	
		  detlogl <- determinant(S,logarithm = TRUE)$modulus[1]
		  kphi <- rep(0,nvar)
		  kphi2 <- rep(0,nlevel)
		  kcount <- matrix(0,nrow=ntint,ncol=nlevel)
		  kmat <- matrix(0,nrow=npoints,ncol=nlevel)

		  foo3 <- .Fortran("ptsamplerwe3",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						ntint	=as.integer(ntint),
						cpar	=as.double(cpar),
						detlogl	=as.double(detlogl),
						z		=as.double(z),
						kphi	=as.integer(kphi),
						kphi2	=as.integer(kphi2),
						kcount  =as.integer(kcount),
						kmat	=as.integer(kmat),
						PACKAGE="DPpackage")	
		
		  kcount <- matrix(foo3$kcount,nrow=ntint,ncol=nlevel)
		  kmat <- matrix(foo3$kmat,nrow=npoints,ncol=nlevel)

		  lgw <- rep(0,npoints)

		  foo4 <- .Fortran("ptsamplerwe4",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						ntint	=as.integer(ntint),
						cpar	=as.double(cpar),
					    z		=as.double(z),
					    detlogl	=as.double(detlogl),
						kcount  =as.integer(kcount),
						kmat	=as.integer(kmat),
					    lgw		=as.double(lgw),
						PACKAGE="DPpackage")	
			  
		  lwg <- foo4$lgw 
		  lwf <- apply(x,1,ltarget) 
		  lr <- lwf - lwg

		#################################
		# Step 4	
		#################################
		  ll1 <- mean(lwg)	  
		  ll2 <- mean(lwf)	  
			  
		  wf <- exp(lwf-ll2)/(sum(exp(lwf-ll2)))	
		  wg <- exp(lwg-ll1)/(sum(exp(lwg-ll1)))
		  L1 <- 0.5*sum(abs(wf-wg))

		#################################
		# Step 5	
		#################################
		  obsmax <- seq(1,npoints)[lr==max(lr)][1]
		  obsmin <- seq(1,npoints)[lr==min(lr)][1]
		
		#################################
		# Step 6	
		#################################

		  ii1 <- rbinom(1,1,tau)	  

	      if(ii1==1)
		  {	  
		     zwork <- z[obsmax,1:nvar]
		     pp1 <- (as.integer((2.0^nlevel)*pnorm(zwork))+runif(nvar))/(2^nlevel)
		     zwork <- qnorm(pp1)	  
		  }
		  if(ii1==0)
		  {
			 zwork <- matrix(rmvnorm(1,mean=z[obsmax,1:nvar],sigma=diag(1,nvar)),nrow=nvar,ncol=1)	  
		  }	  
		  xwork <- m + u%*%zwork
			  
  		  z[obsmin,1:nvar] <- zwork[1:nvar]
		  x[obsmin,1:nvar] <- xwork[1:nvar]

		  kphi <- as.integer((2^nlevel)*pnorm(zwork))
		  kphi2 <- rep(0,nlevel) 
		  for(j1 in 1:nlevel)
		  { 
			  j2 <- nlevel - j1 + 1
			  ind <- sum(2^(j2*(seq(1:nvar)-1))*kphi)
			  kphi2[j2] <- ind + 1
			  kphi <- as.integer(kphi/2)
		  }
		  kmat[obsmin,1:nlevel] <- kphi2[1:nlevel]
		
		#################################
		# Step 7	
		#################################
		  cparc <- rlnorm(1,log(cpar),tune1)	
		
		  qold <- 0
		  qcan <- 0
		
		  foo7.1 <- .Fortran("ptsamplerqe",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						cpar	=as.double(cpar),
						kmat	=as.integer(kmat),
					    eval	=as.double(qold),
						PACKAGE="DPpackage")
		
		  qold <- foo7.1$eval

		  foo7.2 <- .Fortran("ptsamplerqe",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						cpar	=as.double(cparc),
						kmat	=as.integer(kmat),
					    eval	=as.double(qold),
						PACKAGE="DPpackage")
		
		  qcan <- foo7.2$eval
		
		  if(qcan > qold){ cpar <- max(c(minc,cparc))}
		
		#################################
		# Step 8	
		#################################
		  aa <<- list(x=x,L1=L1,cpar=cpar,u=u,uinv=uinv)	
	}


	PTS.sam <- function(aa) function(t)
	{
		#################################
		# Step 0	
		#################################
		  start <- 1
		  end <- nvar
  		  theta <- aa[start:end]

		  start <- end + 1
		  end <- end + nvar*nvar
		  u <- matrix(aa[start:end],nrow=nvar,ncol=nvar)
		
		  start <- end + 1
		  end <- end + nvar*nvar
		  uinv <- matrix(aa[start:end],nrow=nvar,ncol=nvar)
		
		  start <- end + 1
		  end <- end + 1
  		  aratio <- aa[start:end]
		
		  npoints <- nrow(x)
		  nvar <- ncol(x)
		
		  m <- apply(x,2,mean)
		  S <- var(x)
		
		  foo1 <- eigen(S)
		  sqrtLambda <- diag(sqrt(foo1$values))
		  nsqrtLambda <- diag(1/sqrt(foo1$values))
		  M <- foo1$vectors
		  std.points <- function(x)
		  { 
			  uinv%*%(x-m)
		  }	  
		  z <- t(apply(x,1,std.points))

		#################################
		# Step 1	
		#################################
		  A <- matrix(rnorm(nvar*nvar),nrow=nvar,ncol=nvar)
		  out <- qr(A)
		  Q <- qr.Q(out)
		  R <- qr.R(out)
          O <- Q%*%diag(sign(diag(R)))

		  uc <- M%*%sqrtLambda%*%O
	 	  ucinv <- t(M%*%nsqrtLambda%*%O)
			  
		#################################
		# Step 2	
		#################################
		  narea <- 2**nvar
		  massi <- rep(0,narea)
		  mass <- rep(0,narea)
		  parti <- rep(0,nvar)
		  pattern <- rep(0,nvar)
		  patterns <- rep(0,nvar)
		  whicho <- rep(0,npoints)
		  whichn <- rep(0,npoints)	
		  linf <- rep(0,nvar)
		  lsup <- rep(0,nvar)
		  limw <- rep(0,nvar)
		  zwork <- rep(0,nvar)
		  xwork <- rep(0,nvar)

		  foo2 <- .Fortran("ptsamplersam",
						narea	=as.integer(narea),
						nvar	=as.integer(nvar),
						np		=as.integer(npoints),
						nlevel	=as.integer(nlevel),
						cpar	=as.double(cpar),
						m		=as.double(m),
						u		=as.double(uc),
					    z		=as.double(z),
						massi   =as.integer(massi),
						mass	=as.double(mass),
						parti	=as.integer(parti),
						pattern	=as.integer(pattern),
						patterns=as.integer(patterns),
						whicho	=as.integer(whicho),
						whichn	=as.integer(whichn),
						linf	=as.double(linf),
						lsup	=as.double(lsup),
						limw	=as.double(limw),
						zwork	=as.double(zwork),
						xwork	=as.double(xwork),
						PACKAGE="DPpackage")	

		thetac <- foo2$xwork
		thetasc <- foo2$zwork
			
	  #################################
	  # Step 3	
      #################################
		detlogl <- determinant(S,logarithm = TRUE)$modulus[1]
		ntint <- 2**(nvar*nlevel)
		eval <- 0	
		kphi <- rep(0,nvar)
		kphi2 <- rep(0,nlevel)
		kcount <- matrix(0,nrow=ntint,ncol=nlevel)
			
		foo3.1 <- .Fortran("ptsamplermr",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						ntint	=as.integer(ntint),
						cpar	=as.double(cpar),
						detlogl	=as.double(detlogl),
						x		=as.double(x),
						m		=as.double(m),
						uinv	=as.double(uinv),
						kphi	=as.integer(kphi),
						kphi2	=as.integer(kphi2),
						kcount	=as.integer(kcount),
						theta	=as.double(theta),
						eval	=as.double(eval),
						PACKAGE="DPpackage")	
		lpto <- foo3.1$eval	

		foo3.2 <- .Fortran("ptsamplermr",
						np		=as.integer(npoints),
						nvar	=as.integer(nvar),
						nlevel	=as.integer(nlevel),
						ntint	=as.integer(ntint),
						cpar	=as.double(cpar),
						detlogl	=as.double(detlogl),
						x		=as.double(x),
						m		=as.double(m),
						uinv	=as.double(ucinv),
						kphi	=as.integer(kphi),
						kphi2	=as.integer(kphi2),
						kcount	=as.integer(kcount),
						theta	=as.double(thetac),
						eval	=as.double(eval),
						PACKAGE="DPpackage")	
		lptc <- foo3.2$eval	
			
	    lratio <- ltarget(thetac) - ltarget(theta)
		lratio <- lratio + lpto - lptc	

		if(log(runif(1)) < lratio)
		{   
			aratio <- aratio+1
		    theta <- thetac
		    u <- uc
		    uinv <- ucinv
		}   

		aa <<- c(theta,as.vector(u),as.vector(uinv),aratio)
	}

  ######################################################################################
  # call parameters
  ######################################################################################
    cl <- match.call()

  ######################################################################################
  # Support points
  ######################################################################################
	nvar <- dim.theta

    if(status==TRUE)
    {
		if(is.null(support))
		{
			npoints <- 400
			x <- rmvnorm(npoints,mean=rep(0,nvar),sigma=diag(1000,nvar))
		}
		else
		{
			x <- support
			npoints <- nrow(x)
		}
	}
	else
	{
		nvar <- state$dim.theta
		x <- state$support
		npoints <- nrow(x)
	}

  ######################################################################################
  # PTsampler parameters
  ######################################################################################
	
	if(is.null(pts.options))
	{
		nlevel <- 5
		tune1 <- 1
		delta <- 0.2
		max.warmup <- 50000
		minc <- 1
		cpar0 <- 1000
		nadd <- 1000
	}
	else
	{
		nlevel <- pts.options$nlevel
		tune1 <- pts.options$tune1
		delta <- pts.options$delta
		max.warmup <- pts.options$max.warmup
		minc <- pts.options$minc
		cpar0 <- pts.options$cpar0
		nadd <- pts.options$nadd
	}
	tau <- 1

  ######################################################################################
  # MCMC specification
  ######################################################################################

	if(is.null(mcmc))
	{
		nburn <- 1000
		nsave <- 1000
		ndisplay <- 100
	}
	else
	{
		nburn <- mcmc$nburn
		nsave <- mcmc$nsave
		ndisplay <- mcmc$ndisplay
	}

  ######################################################################################
  # Warm-up phase
  ######################################################################################

	if(status==TRUE)
	{
		L1 <- 1
		aa <- list(x=x,L1=L1,cpar=cpar0,u=diag(1,nvar),uinv=diag(1,nvar))
		count <- 0

		cat("\n")
		cat("********************","\n")
		cat("** Warm-up phase  **","\n")
		cat("********************","\n")
		cat("\n")

		while((L1>delta) && (count < max.warmup))
		{
			bb <- sapply(integer(ndisplay),PTS.warmup(aa))
		
			x <- bb[[(ndisplay-1)*5+1]]
			L1 <- bb[[(ndisplay-1)*5+2]]
			cpar <- bb[[(ndisplay-1)*5+3]]
			u <- bb[[(ndisplay-1)*5+4]]
			uinv <- bb[[(ndisplay-1)*5+5]]
			aa <- list(x=x,L1=L1,cpar=cpar,u=u,uinv=uinv)
			count <- count + ndisplay
			cat("Warm-up step #",count,"( L1 = ",round(L1,4),")\n")
		}

		if(L1 <= delta)
		{
			cat("\nConvergence criterion #",L1,"\n")
			cat(nadd," additional steps are now performed\n")
		}
		if(L1 > delta)
		{
			cat("\nConvergence criterion #",L1,"\n")
			cat("Convergence criterion not met\n")
			stop("Try with different initial support points!!\n")
		}

		tau <- 1
		nitr <- nadd
		bb <- sapply(integer(nitr),PTS.warmup(aa))
		x <- bb[[(nitr-1)*5+1]]
		L1 <- bb[[(nitr-1)*5+2]]
		cpar <- bb[[(nitr-1)*5+3]]
		u <- bb[[(nitr-1)*5+4]]
		uinv <- bb[[(nitr-1)*5+5]]

	}
	else
	{
		cpar <- state$cpar
		x <- state$support
		u <- state$u
		uinv <- state$uinv
		L1 <- state$L1
	}

  ######################################################################################
  # Sampling phase
  ######################################################################################
	cat("\n")
    cat("********************","\n")
    cat("** Sampling phase **","\n")
    cat("********************","\n")

  # burn-in

	if(status==FALSE)
	{
		nburn <- 0
		theta <- state$theta
	}
	else
	{
		theta <- x[sample(1:npoints,1),]
	}
	aa <- c(theta,as.vector(u),as.vector(uinv),0)

	if(nburn>0)
	{
		cat("\n")
		count <- 0
		nitr <- min(c(nburn,ndisplay))
		
		while(count < nburn)
		{
			aa <- t(sapply(integer(nitr),PTS.sam(aa), simplify = TRUE))[nitr,]
			count <- count + nitr
			cat("Burn-in step #",count,"\n")
		}
	}	

  # actual MCMC sampling

	thetasave <- NULL
	randsave <- NULL
	acrate <- NULL

    count <- 0
	nitr <- min(c(nsave,ndisplay))
	cat("\n")
	while(count < nsave)
	{
		foo <- t(sapply(integer(nitr),PTS.sam(aa), simplify = TRUE))
		aa <- foo[nitr,]
		count <- count + nitr
		cat("MCMC scan #",count,"\n")
		
		thetasave <- rbind(thetasave,foo[,(1:nvar)])
		start <- nvar+1
		end <- nvar+2*nvar*nvar
		randsave <- rbind(randsave,foo[,(start:end)])
		acrate <- c(acrate,foo[nitr,length(aa)])
	}	

  ######################################################################################
  # Output
  ######################################################################################


	model.name<-"Polya Tree Sampler"		

	acrate <- acrate[length(acrate)]/(nburn+nsave)

	theta <- thetasave[nsave,]
	u <- matrix(randsave[nsave,(1:(nvar*nvar))],nrow=nvar,ncol=nvar)
	uinv <- solve(u)

	state <- list(theta=theta,
				  u=u,
				  uinv=uinv,
				  cpar=cpar,
				  support=x,
				  dim.theta=nvar,
				  L1=L1)

	colnames(thetasave) <- paste("theta",seq(1:nvar),sep="")

	save.state <- list(thetasave=thetasave,
					   randsave=randsave)

	coeff <- apply(thetasave, 2, mean)

	z <- list(call=cl,
			  coefficients=coeff,
			  modelname=model.name,
			  mcmc=mcmc, 
              state=state,
              save.state=save.state,
			  L1=L1,
			  acrate=acrate,
			  dim.theta=dim.theta)
                 
	cat("\n\n")
	class(z) <- c("PTsampler")
	return(z)
}



###                    
### Tools
###
### Copyright: Alejandro Jara, 2010
### Last modification: 21-01-2010.
###


"print.PTsampler" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

	cat("Posterior Inference of Parameters:\n")
	print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)    

	cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.PTsampler" <- function(object, hpd=TRUE, ...) 
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

### 

    dimen1 <- object$dim.theta

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

    names(coef.m) <- names(object$coefficients[1:dimen1])
    names(coef.sd) <- names(object$coefficients[1:dimen1])
    names(coef.se) <- names(object$coefficients[1:dimen1])
    names(coef.l) <- names(object$coefficients[1:dimen1])
    names(coef.u) <- names(object$coefficients[1:dimen1])

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

    ans$acrate <- object$acrate

    class(ans) <- "summaryPTsampler"
    return(ans)
}


"print.summaryPTsampler"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
            
    if (length(x$coefficients)) {
        cat("\nPosterior Inference:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nMH Acceptance Rate = ",x$acrate,"\n")    
    cat("\n\n")
    invisible(x)
}


"plot.PTsampler"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, col="#bdfcc9", ...)
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


   if(is(x, "PTsampler"))
   {
           coef.p <- x$coefficients[1:(x$dim.theta)]
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

}
