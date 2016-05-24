##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
maxim.integrand <- function(y,units.m,mu,Sigma,ID.coords=NULL,poisson.llik) {
	if(length(ID.coords)==0) {
    Sigma.inv <- solve(Sigma)

	integrand <- function(S) {
            if(poisson.llik) {
               llik <- sum(y*S-units.m*exp(S))
            } else {
               llik <- sum(y*S-units.m*log(1+exp(S)))
            }
		diff.S <- S-mu
		q.f_S <- t(diff.S)%*%Sigma.inv%*%(diff.S)
		-0.5*(q.f_S)+
		llik
	}
	
      grad.integrand <- function(S) {
		diff.S <- S-mu
		if(poisson.llik) {
               h <- units.m*exp(S)
            } else {
               h <- units.m*exp(S)/(1+exp(S))
            }
		as.numeric(-Sigma.inv%*%diff.S+(y-h))
	}

	hessian.integrand <- function(S) {
		if(poisson.llik) {
               h1 <- units.m*exp(S)
            } else {
               h1 <- units.m*exp(S)/((1+exp(S))^2)
            }
		res <- -Sigma.inv	
		diag(res) <- diag(res)-h1
		res
	}
	
      out <- list()
	estim <- maxBFGS(function(x) integrand(x),function(x) grad.integrand(x),
	            function(x) hessian.integrand(x),start=mu)
	out$mode <- estim$estimate
	out$Sigma.tilde <- solve(-estim$hessian)             
    
    } else {    	 
    Sigma.inv <- solve(Sigma)
    n.x <- dim(Sigma.inv)[1]
	C.S <- t(sapply(1:n.x,function(i) ID.coords==i))
	integrand <- function(S) {         
		q.f_S <- t(S)%*%Sigma.inv%*%(S)
		eta <- mu+as.numeric(S[ID.coords])
		-0.5*q.f_S+
		sum(y*eta-units.m*log(1+exp(eta)))
	}
	
    grad.integrand <- function(S) {
		eta <- as.numeric(mu+as.numeric(S[ID.coords]))
		h <- units.m*exp(eta)/(1+exp(eta))
		as.numeric(-Sigma.inv%*%S+
		 sapply(1:n.x,function(i) sum((y-h)[C.S[i,]]))	)
	}

	hessian.integrand <- function(S) {
		eta <- as.numeric(mu+as.numeric(S[ID.coords]))
		h <- (units.m*exp(eta)/(1+exp(eta)))
		h1 <- h/(1+exp(eta))
		grad.S.S <-  -Sigma.inv
		diag(grad.S.S) <- diag(grad.S.S)-sapply(1:n.x,function(i) sum(h1[C.S[i,]]))
		as.matrix(grad.S.S)
	}
	out <- list()
	estim <- maxBFGS(integrand,grad.integrand,hessian.integrand,rep(0,n.x))
	out$mode <- estim$estimate
	out$Sigma.tilde <- solve(-estim$hessian)	
    }
    
    return(out)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
maxim.integrand.lr <- function(y,units.m,mu,sigma2,K,poisson.llik) {

	integrand <- function(z) {
            
		s <- as.numeric(K%*%z)
		eta <- as.numeric(mu+s)
            if(poisson.llik) {
               llik <- sum(y*eta-units.m*exp(eta))
            } else {
               llik <- sum(y*eta-units.m*log(1+exp(eta)))
            }
		-0.5*(sum(z^2)/sigma2)+
		llik
	}
	
      grad.integrand <- function(z) {
		s <- as.numeric(K%*%z)
		eta <- as.numeric(mu+s)
            if(poisson.llik) {
               h <- units.m*exp(eta)
            } else {
               h <- units.m*exp(eta)/(1+exp(eta))
            }
		as.numeric(-z/sigma2+t(K)%*%(y-h))		
	}

	hessian.integrand <- function(z) {
		s <- as.numeric(K%*%z)
		eta <- as.numeric(mu+s)
		if(poisson.llik) {
               h1 <- units.m*exp(eta)
            } else {
               h1 <- units.m*exp(eta)/((1+exp(eta))^2)
            }
		out <- -t(K)%*%(K*h1)
		diag(out) <- diag(out)-1/sigma2
		out
	}
	N <- ncol(K)
      out <- list()
	estim <- maxBFGS(function(x) integrand(x),
	                             function(x) grad.integrand(x),
	                             function(x) hessian.integrand(x),rep(0,N))
	out$mode <- estim$estimate
	out$Sigma.tilde <- solve(-estim$hessian)             
            
      return(out)
}


##' @title ID spatial coordinates 
##' @description Creates ID values for the unique set of coordinates.
##' @param data a data frame containing the spatial coordinates. 
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @return a vector of integers indicating the corresponding rows in \code{data} for each distinct coordinate obtained with the \code{\link{unique}} function.
##' @examples
##' x1 <- runif(5)
##' x2 <- runif(5)
##' data <- data.frame(x1=rep(x1,each=3),x2=rep(x2,each=3))
##' ID.coords <- create.ID.coords(data,coords=~x1+x2)
##' data[,c("x1","x2")]==unique(data[,c("x1","x2")])[ID.coords,]
##' 
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export

create.ID.coords <- function(data,coords) {
	if(class(data)!="data.frame") stop("data must be a data frame.")
	if(class(coords)!="formula") stop("coords must a 'formula' object indicating the spatial coordinates in the data.")
    coords <- as.matrix(model.frame(coords,data))
	if(any(is.na(coords))) stop("missing values are not accepted.")
	n <- nrow(coords)
	ID.coords <- rep(NA,n)	
	coords.uni <- unique(coords)
	if(nrow(coords.uni)==n) warning("the unique set of coordinates concides with the provided spatial coordinates in 'coords'.")
	for(i in 1:nrow(coords.uni)) {
		ind <- which(coords.uni[i,1]==coords[,1] & 
		                     coords.uni[i,2]==coords[,2])
        ID.coords[ind] <- i		                     
	}
	return(ID.coords)
}

##' @title Langevin-Hastings MCMC for conditional simulation 
##' @description This function simulates from the conditional distribution of a Gaussian random effect, given binomial or Poisson observations \code{y}.
##' @param mu mean vector of the marginal distribution of the random effect.
##' @param Sigma covariance matrix of the marginal distribution of the random effect.
##' @param y vector of binomial/Poisson observations.
##' @param units.m vector of binomial denominators, or offset if the Poisson model is used.
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param ID.coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID.coords}}. These must be provided if, for example, spatial random effects are defined at household level but some of the covariates are at individual level. \bold{Warning}: the household coordinates must all be distinct otherwise see \code{\link{jitterDupCoords}}. Default is \code{NULL}.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the conditional simulations is displayed. 
##' @param poisson.llik logical; if \code{poisson.llik=TRUE} a Poisson model is used or, if \code{poisson.llik=FALSE}, a binomial model is used.
##' @details \bold{Binomial model.} Conditionally on the random effect \eqn{S}, the data \code{y} follow a binomial distribution with probability \eqn{p} and binomial denominators \code{units.m}. The logistic link function is used for the linear predictor, which assumes the form \deqn{\log(p/(1-p))=S.} 
##' \bold{Poisson model.} Conditionally on the random effect \eqn{S}, the data \code{y} follow a Poisson distribution with mean \eqn{m\lambda}, where \eqn{m} is an offset set through the argument \code{units.m}. The log link function is used for the linear predictor, which assumes the form \deqn{\log(\lambda)=S.} 
##' The random effect \eqn{S} has a multivariate Gaussian distribution with mean \code{mu} and covariance matrix \code{Sigma}. 
##'
##' \bold{Laplace sampling.} This function generates samples from the distribution of \eqn{S} given the data \code{y}. Specifically a Langevin-Hastings algorithm is used to update \eqn{\tilde{S} = \tilde{\Sigma}^{-1/2}(S-\tilde{s})} where \eqn{\tilde{\Sigma}} and \eqn{\tilde{s}} are the inverse of the negative Hessian and the mode of the distribution of \eqn{S} given \code{y}, respectively. At each iteration a new value \eqn{\tilde{s}_{prop}} for \eqn{\tilde{S}} is proposed from a multivariate Gaussian distribution with mean \deqn{\tilde{s}_{curr}+(h/2)\nabla \log f(\tilde{S} | y),}
##' where \eqn{\tilde{s}_{curr}} is the current value for \eqn{\tilde{S}}, \eqn{h} is a tuning parameter and \eqn{\nabla \log f(\tilde{S} | y)} is the the gradient of the log-density of the distribution of \eqn{\tilde{S}} given \code{y}. The tuning parameter \eqn{h} is updated according to the following adaptive scheme: the value of \eqn{h} at the \eqn{i}-th iteration, say \eqn{h_{i}}, is given by \deqn{h_{i} = h_{i-1}+c_{1}i^{-c_{2}}(\alpha_{i}-0.547),}
##' where \eqn{c_{1} > 0} and \eqn{0 < c_{2} < 1} are pre-defined constants, and \eqn{\alpha_{i}} is the acceptance rate at the \eqn{i}-th iteration (\eqn{0.547} is the optimal acceptance rate for a multivariate standard Gaussian distribution). 
##' The starting value for \eqn{h}, and the values for \eqn{c_{1}} and \eqn{c_{2}} can be set through the function \code{\link{control.mcmc.MCML}}.
##' 
##' \bold{Random effects at household-level.} When the data consist of two nested levels, such as households and individuals within households, the argument \code{ID.coords} must be used to define the household IDs for each individual. Let \eqn{i} and \eqn{j} denote the \eqn{i}-th household and the \eqn{j}-th person within that household; the logistic link function then assumes the form \deqn{\log(p_{ij}/(1-p_{ij}))=\mu_{ij}+S_{i}} where the random effects \eqn{S_{i}} are now defined at household level and have mean zero. \bold{Warning:} this modelling option is available only for the binomial model. 
##' @return A list with the following components
##' @return \code{samples}: a matrix, each row of which corresponds to a sample from the predictive distribution.
##' @return \code{h}: vector of the values of the tuning parameter at each iteration of the Langevin-Hastings MCMC algorithm.
##' @seealso \code{\link{control.mcmc.MCML}}, \code{\link{create.ID.coords}}.
##' @examples
##' set.seed(1234)
##' data(data_sim)
##' n.subset <- 50
##' data_subset <- data_sim[sample(1:nrow(data_sim),n.subset),]
##' mu <- rep(0,50)
##' Sigma <- varcov.spatial(coords=data_subset[,c("x1","x2")],
##'               cov.pars=c(1,0.15),kappa=2)$varcov
##' control.mcmc <- control.mcmc.MCML(n.sim=1000,burnin=0,thin=1,
##'                            h=1.65/(n.subset^2/3))               
##' invisible(Laplace.sampling(mu=mu,Sigma=Sigma,
##'                                        y=data_subset$y,units.m=data_subset$units.m,
##'                                        control.mcmc=control.mcmc))
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
Laplace.sampling <- function(mu,Sigma,y,units.m,
                             control.mcmc,ID.coords=NULL,
                             messages=TRUE,
                             plot.correlogram=TRUE,poisson.llik=FALSE) {

   if(poisson.llik==TRUE & length(ID.coords)>0) {
      stop("Two-levels models can only be used for binomial data.")
   }

   if(length(ID.coords)==0) {                                               	
	n.sim <- control.mcmc$n.sim
	n <- length(y)
	S.estim <- maxim.integrand(y,units.m,mu,Sigma,
                                 poisson.llik=poisson.llik)
	Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
	A <- solve(Sigma.sroot)
	Sigma.W.inv <- solve(A%*%Sigma%*%t(A))
	mu.W <- as.numeric(A%*%(mu-S.estim$mode))
	
	cond.dens.W <- function(W,S) {
          if(poisson.llik) {
            llik <- sum(y*S-units.m*exp(S))
          } else {
            llik <- sum(y*S-units.m*log(1+exp(S)))
          }

	    n <- length(y)		
	    diff.W <- W-mu.W
          -0.5*as.numeric(t(diff.W)%*%Sigma.W.inv%*%diff.W)+
          llik
      }

    lang.grad <- function(W,S) {
	    diff.W <- W-mu.W
          if(poisson.llik) {
            h <- as.numeric(units.m*exp(S))
          } else {
            h <- as.numeric(units.m*exp(S)/(1+exp(S)))
          }

	    as.numeric(-Sigma.W.inv%*%diff.W+
          t(Sigma.sroot)%*%(y-h))
    }
  
    h <- control.mcmc$h
    burnin <- control.mcmc$burnin
    thin <- control.mcmc$thin
    c1.h <- control.mcmc$c1.h
    c2.h <- control.mcmc$c2.h
    W.curr <- rep(0,n)
    S.curr <- as.numeric(Sigma.sroot%*%W.curr+S.estim$mode)    
    mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,S.curr))    
    lp.curr <- cond.dens.W(W.curr,S.curr)
    acc <- 0
    sim <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=n)
    if(messages) cat("Conditional simulation (burnin=",control.mcmc$burnin,", thin=",control.mcmc$thin,"): \n",sep="")	
    h.vec <- rep(NA,n.sim)
    for(i in 1:n.sim) {
       W.prop <- mean.curr+h*rnorm(n)
       S.prop <-  as.numeric(Sigma.sroot%*%W.prop+S.estim$mode)
       mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,S.prop))
       lp.prop <- cond.dens.W(W.prop,S.prop)
      
       dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
	   dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
	  
	   log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
	  
	   if(log(runif(1)) < log.prob) {
	      acc <- acc+1
	      W.curr <- W.prop
	      S.curr <- S.prop
		  lp.curr <- lp.prop
		  mean.curr <- mean.prop
	   }

	   if( i > burnin & (i-burnin)%%thin==0) {
	 	 sim[(i-burnin)/thin,] <- S.curr			
	   }
	  
	   h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
       if(messages) cat("Iteration",i,"out of",n.sim,"\r")
       flush.console()
       }
       if(plot.correlogram) {
       acf.plot <- acf(sim[,1],plot=FALSE)
       plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
                 ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")                
       for(i in 2:ncol(sim)) {
          acf.plot <- acf(sim[,i],plot=FALSE)
          lines(acf.plot$lag,acf.plot$acf)
       }
       abline(h=0,lty="dashed",col=2) 
    }
    if(messages) cat("\n")
    } else {
      n.sim <- control.mcmc$n.sim
	  n <- length(y)
	  n.x <- dim(Sigma)[1]
	  S.estim <- maxim.integrand(y,units.m,mu,Sigma,ID.coords)
	  Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
	  A <- solve(Sigma.sroot)
	  Sigma.w.inv <- solve(A%*%Sigma%*%t(A))
	  mu.w <- -as.numeric(A%*%S.estim$mode)	
	  
	  cond.dens.W <- function(W,S) {
	    diff.w <- W-mu.w
	    eta <- mu+as.numeric(S[ID.coords])		
        -0.5*as.numeric(t(diff.w)%*%Sigma.w.inv%*%diff.w)+
        sum(y*eta-units.m*log(1+exp(eta)))
      }

      lang.grad <- function(W,S) {
	    diff.w <- W-mu.w
	    eta <- mu+as.numeric(S[ID.coords])
        der <- units.m*exp(eta)/(1+exp(eta))     
        grad.S <- sapply(1:n.x,function(i) sum((y-der)[C.S[i,]]))
	    as.numeric(-Sigma.w.inv%*%(W-mu.w)+
        t(Sigma.sroot)%*%grad.S)
      }
      
      h <- control.mcmc$h
      burnin <- control.mcmc$burnin
      thin <- control.mcmc$thin
      c1.h <- control.mcmc$c1.h
      c2.h <- control.mcmc$c2.h
      W.curr <- rep(0,n.x)
      S.curr <- as.numeric(Sigma.sroot%*%W.curr+S.estim$mode) 
      C.S <- t(sapply(1:n.x,function(i) ID.coords==i))   
      mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,S.curr))    
      lp.curr <- cond.dens.W(W.curr,S.curr)
      acc <- 0
      sim <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=n.x)
      if(messages) cat("Conditional simulation (burnin=",
      control.mcmc$burnin,", thin=",control.mcmc$thin,"): \n",sep="")	
      h.vec <- rep(NA,n.sim)
      for(i in 1:n.sim) {
       W.prop <- mean.curr+h*rnorm(n.x)
       S.prop <-  as.numeric(Sigma.sroot%*%W.prop+S.estim$mode)
       mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,S.prop))
       lp.prop <- cond.dens.W(W.prop,S.prop)
      
       dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
	   dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
	  
	   log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
	  
	   if(log(runif(1)) < log.prob) {
	      acc <- acc+1
	      W.curr <- W.prop
	      S.curr <- S.prop
		lp.curr <- lp.prop
		mean.curr <- mean.prop
	   }

	   if( i > burnin & (i-burnin)%%thin==0) {
	 	 sim[(i-burnin)/thin,] <- S.curr			
	   }
	  
	   h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
       if(messages) cat("Iteration",i,"out of",n.sim,"\r")
       flush.console()
       }
        if(plot.correlogram) {
       acf.plot <- acf(sim[,1],plot=FALSE)
       plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
                 ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")                
       for(i in 2:ncol(sim)) {
          acf.plot <- acf(sim[,i],plot=FALSE)
          lines(acf.plot$lag,acf.plot$acf)
       }
       abline(h=0,lty="dashed",col=2)           
       }
       if(messages) cat("\n")
    }
    out.sim <- list(samples=sim,h=h.vec)
    class(out.sim) <- "mcmc.PrevMap"
    return(out.sim)           
}

##' @title Langevin-Hastings MCMC for conditional simulation (low-rank approximation)
##' @description This function simulates from the conditional distribution of the random effects of binomial and Poisson models.
##' @param mu mean vector of the linear predictor.
##' @param sigma2 variance of the random effect.
##' @param K random effect design matrix, or kernel matrix for the low-rank approximation.
##' @param y vector of binomial/Poisson observations.
##' @param units.m vector of binomial denominators, or offset if the Poisson model is used.
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the conditional simulations is displayed.
##' @param poisson.llik logical; if \code{poisson.llik=TRUE} a Poisson model is used or, if \code{poisson.llik=FALSE}, a binomial model is used.
##' @details \bold{Binomial model.} Conditionally on \eqn{Z}, the data \code{y} follow a binomial distribution with probability \eqn{p} and binomial denominators \code{units.m}. Let \eqn{K} denote the random effects design matrix; a logistic link function is used, thus the linear predictor assumes the form \deqn{\log(p/(1-p))=\mu + KZ} where \eqn{\mu} is the mean vector component defined through \code{mu}.
##' \bold{Poisson model.} Conditionally on \eqn{Z}, the data \code{y} follow a Poisson distribution with mean \eqn{m\lambda}, where \eqn{m} is an offset set through the argument \code{units.m}. Let \eqn{K} denote the random effects design matrix; a log link function is used, thus the linear predictor assumes the form \deqn{\log(\lambda)=\mu + KZ} where \eqn{\mu} is the mean vector component defined through \code{mu}.
##' The random effect \eqn{Z} has iid components distributed as zero-mean Gaussian variables with variance \code{sigma2}. 
##'
##' \bold{Laplace sampling.} This function generates samples from the distribution of \eqn{Z} given the data \code{y}. Specifically, a Langevin-Hastings algorithm is used to update \eqn{\tilde{Z} = \tilde{\Sigma}^{-1/2}(Z-\tilde{z})} where \eqn{\tilde{\Sigma}} and \eqn{\tilde{z}} are the inverse of the negative Hessian and the mode of the distribution of \eqn{Z} given \code{y}, respectively. At each iteration a new value \eqn{\tilde{z}_{prop}} for \eqn{\tilde{Z}} is proposed from a multivariate Gaussian distribution with mean \deqn{\tilde{z}_{curr}+(h/2)\nabla \log f(\tilde{Z} | y),}
##' where \eqn{\tilde{z}_{curr}} is the current value for \eqn{\tilde{Z}}, \eqn{h} is a tuning parameter and \eqn{\nabla \log f(\tilde{Z} | y)} is the the gradient of the log-density of the distribution of \eqn{\tilde{Z}} given \code{y}. The tuning parameter \eqn{h} is updated according to the following adaptive scheme: the value of \eqn{h} at the \eqn{i}-th iteration, say \eqn{h_{i}}, is given by \deqn{h_{i} = h_{i-1}+c_{1}i^{-c_{2}}(\alpha_{i}-0.547),}
##' where \eqn{c_{1} > 0} and \eqn{0 < c_{2} < 1} are pre-defined constants, and \eqn{\alpha_{i}} is the acceptance rate at the \eqn{i}-th iteration (\eqn{0.547} is the optimal acceptance rate for a multivariate standard Gaussian distribution). 
##' The starting value for \eqn{h}, and the values for \eqn{c_{1}} and \eqn{c_{2}} can be set through the function \code{\link{control.mcmc.MCML}}.
##' 
##' @return A list with the following components
##' @return \code{samples}: a matrix, each row of which corresponds to a sample from the predictive distribution.
##' @return \code{h}: vector of the values of the tuning parameter at each iteration of the Langevin-Hastings MCMC algorithm.
##' @seealso \code{\link{control.mcmc.MCML}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
Laplace.sampling.lr <- function(mu,sigma2,K,y,units.m,
                                control.mcmc,
                                messages=TRUE,
                                plot.correlogram=TRUE,
                                poisson.llik=FALSE) {
	n.sim <- control.mcmc$n.sim
	n <- length(y)
	N <- ncol(K)
	S.estim <- maxim.integrand.lr(y,units.m,mu,sigma2,K,
                 poisson.llik=poisson.llik)
	Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
	A <- solve(Sigma.sroot)
	Sigma.W.inv <- solve(sigma2*A%*%t(A))
	mu.W <- -as.numeric(A%*%S.estim$mode)
	
      cond.dens.W <- function(W,Z) {
	   diff.w <- W-mu.W
	   S <- as.numeric(K%*%Z)
         eta <- as.numeric(mu+S)
         if(poisson.llik) {
            llik <- sum(y*eta-units.m*exp(eta))
         } else {
            llik <- sum(y*eta-units.m*log(1+exp(eta)))
         }

         -0.5*as.numeric(t(diff.w)%*%Sigma.W.inv%*%diff.w)+
         llik
      }

      lang.grad <- function(W,Z) {
 	   diff.w <- W-mu.W
 	   S <- as.numeric(K%*%Z)
         eta <- as.numeric(mu+S)
         if(poisson.llik) {
            h <- units.m*exp(eta)
         } else {
            h <- units.m*exp(eta)/(1+exp(eta)) 
         }

         grad.z <- t(K)%*%(y-h)
         as.numeric(-Sigma.W.inv%*%(W-mu.W)+
         t(Sigma.sroot)%*%c(grad.z))
      }

    h <- control.mcmc$h    
    burnin <- control.mcmc$burnin
    thin <- control.mcmc$thin
    c1.h <- control.mcmc$c1.h
    c2.h <- control.mcmc$c2.h
    W.curr <- rep(0,N)
    Z.curr <- as.numeric(Sigma.sroot%*%W.curr+S.estim$mode)    
    mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,Z.curr))    
    lp.curr <- cond.dens.W(W.curr,Z.curr)
    acc <- 0
    sim <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=N)
    if(messages) cat("Conditional simulation (burnin=",control.mcmc$burnin,", thin=",control.mcmc$thin,"): \n",sep="")	
    h.vec <- rep(NA,n.sim)
    for(i in 1:n.sim) {
       W.prop <- mean.curr+h*rnorm(N)
       Z.prop <-  as.numeric(Sigma.sroot%*%W.prop+S.estim$mode)
       mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,Z.prop))
       lp.prop <- cond.dens.W(W.prop,Z.prop)
      
       dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
	   dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
	  
	   log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
	  
	   if(log(runif(1)) < log.prob) {
	      acc <- acc+1
	      W.curr <- W.prop
	      Z.curr <- Z.prop
		  lp.curr <- lp.prop
		  mean.curr <- mean.prop
	   }

	   if( i > burnin & (i-burnin)%%thin==0) {
	 	 sim[(i-burnin)/thin,] <- Z.curr			
	   }
	  
	   h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
       if(messages) cat("Iteration",i,"out of",n.sim,"\r")
       flush.console()
       }
       if(plot.correlogram) {
          acf.plot <- acf(sim[,1],plot=FALSE)
          plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
                 ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")                
          for(i in 2:ncol(sim)) {
          	 acf.plot <- acf(sim[,i],plot=FALSE)
             lines(acf.plot$lag,acf.plot$acf)
          }
          abline(h=0,lty="dashed",col=2) 
       }
       if(messages) cat("\n")       
       out.sim <- list(samples=sim,h=h.vec)
       class(out.sim) <- "mcmc.PrevMap"
       return(out.sim)        
}

##' @title Control settings for the MCMC algorithm used for classical inference on a binomial logistic model
##' @description This function defines the options for the MCMC algorithm used in the Monte Carlo maximum likelihood method.
##' @param n.sim number of simulations.
##' @param burnin length of the burn-in period. 
##' @param thin only every \code{thin} iterations, a sample is stored; default is \code{thin=1}.
##' @param h tuning parameter of the proposal distribution; default is \code{h=0.05}.
##' @param c1.h value of \eqn{c_{1}} used in the adaptive scheme for \code{h}; default is \code{c1.h=0.01}. See also 'Details' in \code{\link{binomial.logistic.MCML}}
##' @param c2.h value of \eqn{c_{2}} used in the adaptive scheme for \code{h}; default is \code{c1.h=0.01}. See also 'Details' in \code{\link{binomial.logistic.MCML}}
##' @return A list with processed arguments to be passed to the main function.
##' @examples
##' control.mcmc <- control.mcmc.MCML(n.sim=1000,burnin=100,thin=1,h=0.05)
##' str(control.mcmc)
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
control.mcmc.MCML <- function(n.sim,burnin,thin=1,h,c1.h=0.01,c2.h=0.0001) {
	if(n.sim < burnin) stop("n.sim cannot be smaller than burnin.")
	if(thin <= 0) stop("thin must be positive")
	if((n.sim-burnin)%%thin!=0) stop("thin must be a divisor of (n.sim-burnin)")
	if(h < 0) stop("h must be positive.")
	if(c1.h < 0) stop("c1.h must be positive.")
	if(c2.h < 0 | c2.h > 1) stop("c2.h must be between 0 and 1.")
	res <- list(n.sim=n.sim,burnin=burnin,thin=thin,h=h,c1.h=c1.h,c2.h=c2.h)
	class(res) <- "mcmc.MCML.PrevMap"
	return(res)
}

##' @title Matern kernel
##' @description This function computes values of the Matern kernel for given distances and parameters.
##' @param u a vector, matrix or array with values of the distances between pairs of data locations.
##' @param rho value of the (re-parametrized) scale parameter; this corresponds to the re-parametrization \code{rho = 2*sqrt(kappa)*phi}. 
##' @param kappa value of the shape parameter.
##' @details The Matern kernel is defined as:
##' \deqn{
##' K(u; \phi, \kappa) = \frac{\Gamma(\kappa + 1)^{1/2}\kappa^{(\kappa+1)/4}u^{(\kappa-1)/2}}{\pi^{1/2}\Gamma((\kappa+1)/2)\Gamma(\kappa)^{1/2}(2\kappa^{1/2}\phi)^{(\kappa+1)/2}}\mathcal{K}_{\kappa}(u/\phi), u > 0,	
##' }	
##' where \eqn{\phi} and \eqn{\kappa} are the scale and shape parameters, respectively, and \eqn{\mathcal{K}_{\kappa}(.)} is the modified Bessel function of the third kind of order \eqn{\kappa}. The family is valid for \eqn{\phi > 0} and \eqn{\kappa > 0}. 
##' @return  A vector matrix or array, according to the argument u, with the values of the Matern kernel function for the given distances.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export 
matern.kernel <- function(u,rho,kappa) {
	u <- u+1e-100
	if(kappa==2) {
		out <- 4*exp(-2*sqrt(2)*u/rho)/((pi^0.5)*rho)
	} else  {
		out <- (2*gamma(kappa+1)^(0.5))*((kappa)^((kappa+1)/4))*
		          (u^((kappa-1)/2))/((pi^0.5)*gamma((kappa+1)/2)*
		           (gamma(kappa)^0.5)*(rho^((kappa+1)/2)))*
		           besselK(2*sqrt(kappa)*u/rho,(kappa-1)/2)
	}
    out		
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom geoR varcov.spatial
##' @importFrom maxLik maxBFGS
geo.MCML <- function(formula,units.m,coords,data,ID.coords,
                                                   par0,control.mcmc,kappa,
                                                   fixed.rel.nugget,
                                                   start.cov.pars,
                                                   method,messages,
                                                   plot.correlogram,
                                                   poisson.llik) {
    start.cov.pars <- as.numeric(start.cov.pars)
    if(any(start.cov.pars < 0)) stop("start.cov.pars must be positive")
    if((length(fixed.rel.nugget)==0 & length(start.cov.pars)!=2) |
       (length(fixed.rel.nugget)>0 & length(start.cov.pars)!=1)) stop("wrong values for start.cov.pars")
    kappa <- as.numeric(kappa)
    if(any(is.na(data))) stop("missing values are not accepted")     
        
    der.phi <- function(u,phi,kappa) {
        u <- u+10e-16
        if(kappa==0.5) {
           out <- (u*exp(-u/phi))/phi^2
        } else {
           out <- ((besselK(u/phi,kappa+1)+besselK(u/phi,kappa-1))*
                  phi^(-kappa-2)*u^(kappa+1))/(2^kappa*gamma(kappa))-
                  (kappa*2^(1-kappa)*besselK(u/phi,kappa)*phi^(-kappa-1)*
                  u^kappa)/gamma(kappa)
        }
        out
     }

    der2.phi <- function(u,phi,kappa) {
        u <- u+10e-16
        if(kappa==0.5) {
        	out <- (u*(u-2*phi)*exp(-u/phi))/phi^4
        } else {
            bk <- besselK(u/phi,kappa)
            bk.p1 <- besselK(u/phi,kappa+1)
            bk.p2 <- besselK(u/phi,kappa+2)
            bk.m1 <- besselK(u/phi,kappa-1)
            bk.m2 <- besselK(u/phi,kappa-2)
            out <- (2^(-kappa-1)*phi^(-kappa-4)*u^kappa*(bk.p2*u^2+2*bk*u^2+
                bk.m2*u^2-4*kappa*bk.p1*phi*u-4*
                bk.p1*phi*u-4*kappa*bk.m1*phi*u-4*bk.m1*phi*u+
                4*kappa^2*bk*phi^2+4*kappa*bk*phi^2))/(gamma(kappa))
        }
        out
    }
    
    matern.grad.phi <- function(U,phi,kappa) {
    	n <- attr(U,"Size")
        grad.phi.mat <- matrix(NA,nrow=n,ncol=n)
        ind <- lower.tri(grad.phi.mat)             
        grad.phi <- der.phi(as.numeric(U),phi,kappa)
        grad.phi.mat[ind] <-  grad.phi
        grad.phi.mat <- t(grad.phi.mat)
        grad.phi.mat[ind] <-  grad.phi
        diag(grad.phi.mat) <- rep(der.phi(0,phi,kappa),n)
        grad.phi.mat
     }

     matern.hessian.phi <- function(U,phi,kappa) {
     	n <- attr(U,"Size")
	    hess.phi.mat <- matrix(NA,nrow=n,ncol=n)
        ind <- lower.tri(hess.phi.mat)             
        hess.phi <- der2.phi(as.numeric(U),phi,kappa)
        hess.phi.mat[ind] <-  hess.phi
        hess.phi.mat <- t(hess.phi.mat)
        hess.phi.mat[ind] <-  hess.phi
        diag(hess.phi.mat) <- rep(der2.phi(0,phi,kappa),n)
        hess.phi.mat
     }
     
    if(length(ID.coords)==0) {	
    	coords <- as.matrix(model.frame(coords,data))
	if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("wrong set of coordinates.")

    mf <- model.frame(formula,data=data)
    y <- as.numeric(model.response(mf))
    n <- length(y)
    if(poisson.llik && length(units.m)==0) {
      units.m <- rep(1,n)
    } else {
      units.m <-  as.numeric(model.frame(units.m,data)[,1])  
    }
     
    if(is.numeric(units.m)==FALSE) stop("'units.m' must be numeric.")
    
    D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
    p <- ncol(D)
    if(any(par0[-(1:p)] <= 0)) stop("the covariance parameters in 'par0' must be positive.")
    beta0 <- par0[1:p]
    mu0 <- as.numeric(D%*%beta0)
		
    sigma2.0 <- par0[p+1]
    phi0 <- par0[p+2]
    if(length(fixed.rel.nugget)>0){
	   if(length(fixed.rel.nugget) != 1 | fixed.rel.nugget < 0) stop("negative fixed nugget value or wrong length")
	   if(length(par0)!=(p+2)) stop("wrong length of par0")	
	   tau2.0 <- fixed.rel.nugget	
	   cat("Fixed relative variance of the nugget effect:",tau2.0,"\n")	
	} else {
	   if(length(par0)!=(p+3)) stop("wrong length of par0")		
	   tau2.0 <- par0[p+3]	
	}	
	
	U <- dist(coords)
	Sigma0 <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(sigma2.0,phi0),
	                                        nugget=tau2.0,kappa=kappa)$varcov		

    n.sim <- control.mcmc$n.sim
    S.sim.res <- Laplace.sampling(mu0,Sigma0,y,units.m,control.mcmc,
	             plot.correlogram=plot.correlogram,messages=messages,
                   poisson.llik=poisson.llik)
    S.sim <- S.sim.res$samples
	
    log.integrand <- function(S,val) {
       if(poisson.llik) {
          llik <- sum(y*S-units.m*exp(S))
       } else {
          llik <- sum(y*S-units.m*log(1+exp(S)))
       }
       diff.S <- S-val$mu
       q.f_S <-    t(diff.S)%*%val$R.inv%*%(diff.S)
       -0.5*(n*log(val$sigma2)+
       val$ldetR+
       q.f_S/val$sigma2)+
       llik
    }	
   
    compute.log.f <- function(par,ldetR=NA,R.inv=NA) {            
       beta <- par[1:p]       
       phi <- exp(par[p+2])
       if(length(fixed.rel.nugget)>0) {
       	  nu2 <- fixed.rel.nugget
       } else {
          nu2 <- exp(par[p+3])
       }
       val <- list()
       val$mu <- as.numeric(D%*%beta)       
       val$sigma2 <- exp(par[p+1])
       if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {                              
   	      R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov
          val$ldetR <- determinant(R)$modulus
          val$R.inv <- solve(R)
       } else {
   	      val$ldetR <- ldetR
          val$R.inv <- R.inv 
       }    
       sapply(1:(dim(S.sim)[1]),function(i) log.integrand(S.sim[i,],val))
     }      

     log.f.tilde <- compute.log.f(c(beta0,log(c(sigma2.0,phi0,tau2.0/sigma2.0))))

     MC.log.lik <- function(par) {
        log(mean(exp(compute.log.f(par)-log.f.tilde)))
     }          
    
     grad.MC.log.lik <- function(par) {
        beta <- par[1:p]; mu <- D%*%beta
        sigma2 <- exp(par[p+1])
        phi <- exp(par[p+2])
        if(length(fixed.rel.nugget) > 0) {
          nu2 <- fixed.rel.nugget	
        } else {
        	  nu2 <- exp(par[p+3])	
        }        

        R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov
        R.inv <- solve(R)
        ldetR <- determinant(R)$modulus                            
                             
        exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)
        L.m <- sum(exp.fact)
        exp.fact <- exp.fact/L.m           
        
        R1.phi <- matern.grad.phi(U,phi,kappa)                        
        m1.phi <- R.inv%*%R1.phi
        t1.phi <- -0.5*sum(diag(m1.phi))
        m2.phi <- m1.phi%*%R.inv; rm(m1.phi)

        if(length(fixed.rel.nugget)==0) {
        t1.nu2 <- -0.5*sum(diag(R.inv))
        m2.nu2 <- R.inv%*%R.inv
        }
                                    
        gradient.S <- function(S) {
      
           diff.S <- S-mu
           q.f <- t(diff.S)%*%R.inv%*%diff.S
           grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2

           grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2           

           grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.S)%*%m2.phi%*%(diff.S))/sigma2)*phi
           
           if(length(fixed.rel.nugget)==0){
           	 grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.S)%*%m2.nu2%*%(diff.S))/sigma2)*nu2
             out <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2)
           } else {
           	 out <- c(grad.beta,grad.log.sigma2,grad.log.phi)
           }        
           return(out)   
        }
        out <- rep(0,length(par))
        for(i in 1:(dim(S.sim)[1])) {
   	        out <- out + exp.fact[i]*gradient.S(S.sim[i,])
        }
        out
      }
      
      hess.MC.log.lik <- function(par) {
         beta <- par[1:p]; mu <- D%*%beta
         sigma2 <- exp(par[p+1])
         phi <- exp(par[p+2])
        if(length(fixed.rel.nugget) > 0) {
          nu2 <- fixed.rel.nugget	
        } else {
        	  nu2 <- exp(par[p+3])	
        }     

         R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov
         R.inv <- solve(R)
         ldetR <- determinant(R)$modulus                            
                             
         exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)
         L.m <- sum(exp.fact)
         exp.fact <- exp.fact/L.m           
        
         R1.phi <- matern.grad.phi(U,phi,kappa)                        
         m1.phi <- R.inv%*%R1.phi
         t1.phi <- -0.5*sum(diag(m1.phi))
         m2.phi <- m1.phi%*%R.inv; rm(m1.phi)
     
         if(length(fixed.rel.nugget)==0) {    
         t1.nu2 <- -0.5*sum(diag(R.inv))
         m2.nu2 <- R.inv%*%R.inv   
         t2.nu2 <- 0.5*sum(diag(m2.nu2))
         n2.nu2 <- 2*R.inv%*%m2.nu2
         t2.nu2.phi <- 0.5*sum(diag(R.inv%*%R1.phi%*%R.inv))
         n2.nu2.phi <- R.inv%*%(R.inv%*%R1.phi+
                               R1.phi%*%R.inv)%*%R.inv  
         }
         
         R2.phi <- matern.hessian.phi(U,phi,kappa)
         t2.phi <- -0.5*sum(diag(R.inv%*%R2.phi-R.inv%*%R1.phi%*%R.inv%*%R1.phi))
         n2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv          
         
         H <- matrix(0,nrow=length(par),ncol=length(par)) 
         H[1:p,1:p] <- -t(D)%*%R.inv%*%D/sigma2                       
         hessian.S <- function(S,ef) {
    
           diff.S <- S-mu
           q.f <- t(diff.S)%*%R.inv%*%diff.S
           grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2

           grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
           grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.S)%*%m2.phi%*%(diff.S))/sigma2)*phi
           if(length(fixed.rel.nugget)==0) {
           grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.S)%*%m2.nu2%*%(diff.S))/sigma2)*nu2           
           g <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2)	
           } else {
           g <- c(grad.beta,grad.log.sigma2,grad.log.phi)		
           }
                 
           grad2.log.beta.lsigma2 <- -t(D)%*%R.inv%*%(diff.S)/sigma2
     
           grad2.log.beta.lphi <- -phi*as.numeric(t(D)%*%m2.phi%*%(diff.S))/sigma2
         
           
     
           grad2.log.lsigma2.lsigma2 <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
                                                    grad.log.sigma2
   
           grad2.log.lphi.lphi <-(t2.phi-0.5*t(diff.S)%*%n2.phi%*%(diff.S)/sigma2)*phi^2+
                                            grad.log.phi         
           
           if(length(fixed.rel.nugget)==0) {
               grad2.log.beta.lnu2 <- -nu2*as.numeric(t(D)%*%m2.nu2%*%(diff.S))/sigma2           
               grad2.log.lnu2.lnu2 <- (t2.nu2-0.5*t(diff.S)%*%n2.nu2%*%(diff.S)/sigma2)*nu2^2+
                                         grad.log.nu2
               grad2.log.lnu2.lphi <- (t2.nu2.phi-0.5*t(diff.S)%*%n2.nu2.phi%*%(diff.S)/sigma2)*phi*nu2   
               H[1:p,p+1] <- H[p+1,1:p]<- grad2.log.beta.lsigma2
               H[1:p,p+2] <- H[p+2,1:p]<- grad2.log.beta.lphi
               H[1:p,p+3] <- H[p+3,1:p]<- grad2.log.beta.lnu2
   	           H[p+1,p+1] <-  grad2.log.lsigma2.lsigma2
   	           H[p+1,p+3] <- H[p+3,p+1] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
   	           H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
   	           H[p+3,p+3] <- grad2.log.lnu2.lnu2
   	           H[p+2,p+3] <- H[p+3,p+2] <- grad2.log.lnu2.lphi
   	           H[p+2,p+2] <- grad2.log.lphi.lphi
           } else {
           	   H[1:p,p+1] <- H[p+1,1:p]<- grad2.log.beta.lsigma2
               H[1:p,p+2] <- H[p+2,1:p]<- grad2.log.beta.lphi
   	           H[p+1,p+1] <-  grad2.log.lsigma2.lsigma2
   	           H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
   	           H[p+2,p+2] <- grad2.log.lphi.lphi
           }                                            
     
   	       out <- list()
   	       out$mat1<- ef*(g%*%t(g)+H)                             
   	       out$g <- g*ef
   	       out
         }            
             
         a <- rep(0,length(par))
         A <- matrix(0,length(par),length(par))
        for(i in 1:(dim(S.sim)[1])) {
   	       out.i <- hessian.S(S.sim[i,],exp.fact[i])
   	       a <- a+out.i$g
           A <- A+out.i$mat1
        }  
       (A-a%*%t(a))
      }
      
      if(messages) cat("Estimation: \n")
      start.par <- c(par0[1:(p+1)],start.cov.pars)
      start.par[-(1:p)] <- log(start.par[-(1:p)])
      
      estim <- list()
      if(method=="BFGS") {
         estimBFGS <- maxBFGS(MC.log.lik,grad.MC.log.lik,hess.MC.log.lik,
                           start.par,print.level=1*messages)   
         estim$estimate <- estimBFGS$estimate  
         estim$covariance <- solve(-estimBFGS$hessian)                
         estim$log.lik <- estimBFGS$maximum                                                 
      }
      
      if(method=="nlminb") {
      	estimNLMINB <- nlminb(start.par,function(x) -MC.log.lik(x),
      	            function(x) -grad.MC.log.lik(x),
      	            function(x) -hess.MC.log.lik(x),control=list(trace=1*messages)) 
      	estim$estimate <- estimNLMINB$par            
      	estim$covariance <- solve(-hess.MC.log.lik(estimNLMINB$par))
      	estim$log.lik <- -estimNLMINB$objective              
      }
      
      } else {  
      coords <- unique(as.matrix(model.frame(coords,data)))
	  if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("wrong set of coordinates.")
    	
      mf <- model.frame(formula,data=data)
	  y <- as.numeric(model.response(mf))
	  n <- length(y)
	  n.x <- nrow(coords)
	  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
	  p <- ncol(D)
	  beta0 <- par0[1:p]
	  mu0 <- as.numeric(D%*%beta0)

	sigma2.0 <- par0[p+1]
	phi0 <- par0[p+2]
    if(length(fixed.rel.nugget)>0){
	   if(length(fixed.rel.nugget) != 1 | fixed.rel.nugget < 0) stop("negative fixed nugget value or wrong length")
	   if(length(par0)!=(p+2)) stop("wrong length of par0")	
	   tau2.0 <- fixed.rel.nugget
	   cat("Fixed relative variance of the nugget effect:",tau2.0,"\n")	
	} else {
	   if(length(par0)!=(p+3)) stop("wrong length of par0")		
	   tau2.0 <- par0[p+3]		   
	}	
	U <- dist(coords)
	Sigma0 <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(sigma2.0,phi0),
	                                        nugget=tau2.0,kappa=kappa)$varcov		
	                                        	
    n.sim <- control.mcmc$n.sim
	S.sim.res <- Laplace.sampling(mu0,Sigma0,y,units.m,control.mcmc,ID.coords,
	               plot.correlogram=plot.correlogram,messages=messages,
                     poisson.llik=FALSE)
    S.sim <- S.sim.res$samples
    log.integrand <- function(S,val) {
        n <- length(y)
      n.x <- length(S)
      eta <- S[ID.coords]+val$mu 
      q.f_S <- t(S)%*%val$R.inv%*%S/val$sigma2   
      -0.5*(n.x*log(val$sigma2)+val$ldetR+q.f_S)+
      sum(y*eta-units.m*log(1+exp(eta)))
   }

   compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    if(length(fixed.rel.nugget)>0) {
       	  nu2 <- fixed.rel.nugget
    } else {
          nu2 <- exp(par[p+3])
    }
    phi <- exp(par[p+2])
    val <- list()  
    val$sigma2 <- sigma2
    val$mu <- as.numeric(D%*%beta)   
    if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {                              
   	  R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov                       
      val$ldetR <- determinant(R)$modulus
      val$R.inv <- solve(R)
     } else {
   	  val$ldetR <- ldetR
      val$R.inv <- R.inv 
    }    
      sapply(1:(dim(S.sim)[1]),function(i) log.integrand(S.sim[i,],val))
   } 
   
   log.f.tilde <- compute.log.f(c(beta0,log(c(sigma2.0,phi0,tau2.0/sigma2.0))))
     
   MC.log.lik <- function(par) {
        log(mean(exp(compute.log.f(par)-log.f.tilde)))
   }

   grad.MC.log.lik <- function(par) {      
      beta <- par[1:p]; mu <- D%*%beta
      sigma2 <- exp(par[p+1])      
      phi <- exp(par[p+2])
      
      if(length(fixed.rel.nugget)==0) {
      	nu2 <- exp(par[p+3])
      } else {
      	nu2 <- fixed.rel.nugget
      }
   
      R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov                              
      R.inv <- solve(R)
      ldetR <- determinant(R)$modulus                            
                             
      exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)    
      L.m <- sum(exp.fact)
      exp.fact <- exp.fact/L.m          

      R1.phi <- matern.grad.phi(U,phi,kappa)                           
      m1.phi <- R.inv%*%R1.phi
      t1.phi <- -0.5*sum(diag(m1.phi))
      m2.phi <- m1.phi%*%R.inv; rm(m1.phi)
      
      if(length(fixed.rel.nugget)==0){                          
      t1.nu2 <- -0.5*sum(diag(R.inv))
      m2.nu2 <- R.inv%*%R.inv
      }
                                      
      gradient.S <- function(S) {
          eta <- mu+S[ID.coords]
          h <- units.m*exp(eta)/(1+exp(eta))
          q.f <- t(S)%*%R.inv%*%S
      
          grad.beta <-  t(D)%*%(y-h)

          grad.log.sigma2 <- (-n.x/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2

          grad.log.phi <- (t1.phi+0.5*as.numeric(t(S)%*%m2.phi%*%(S))/sigma2)*phi
          
          if(length(fixed.rel.nugget)==0) {
            grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(S)%*%m2.nu2%*%(S))/sigma2)*nu2          	
          	out <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2) 
          } else {
          	out <- c(grad.beta,grad.log.sigma2,grad.log.phi) 
          }
          out                    
       }
       out <- rep(0,length(par))
       for(i in 1:(dim(S.sim)[1])) {
   	      out <- out + exp.fact[i]*gradient.S(S.sim[i,])
       }
       out
   }
        
   hess.MC.log.lik <- function(par) {
      beta <- par[1:p]; mu <- D%*%beta
      sigma2 <- exp(par[p+1])
      if(length(fixed.rel.nugget)==0) {
      	nu2 <- exp(par[p+3])
      } else {
      	nu2 <- fixed.rel.nugget
      }
      phi <- exp(par[p+2])
   
      R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                     cov.pars=c(1,phi),
	                                     nugget=nu2,kappa=kappa)$varcov                              
      R.inv <- solve(R)
      ldetR <- determinant(R)$modulus                            
                             
      exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)    
      L.m <- sum(exp.fact)
      exp.fact <- exp.fact/L.m          

      R1.phi <- matern.grad.phi(U,phi,kappa)                           
      m1.phi <- R.inv%*%R1.phi
      t1.phi <- -0.5*sum(diag(m1.phi))
      m2.phi <- m1.phi%*%R.inv; rm(m1.phi)
      
      if(length(fixed.rel.nugget)==0){                          
      t1.nu2 <- -0.5*sum(diag(R.inv))
      m2.nu2 <- R.inv%*%R.inv         
      t2.nu2 <- 0.5*sum(diag(m2.nu2))
      n2.nu2 <- 2*R.inv%*%m2.nu2      
      t2.nu2.phi <- 0.5*sum(diag(R.inv%*%R1.phi%*%R.inv))
      n2.nu2.phi <- R.inv%*%(R.inv%*%R1.phi+
                            R1.phi%*%R.inv)%*%R.inv 
      }

      R2.phi <- matern.hessian.phi(U,phi,kappa)
      t2.phi <- -0.5*sum(diag(R.inv%*%R2.phi-R.inv%*%R1.phi%*%R.inv%*%R1.phi))
      n2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv          
   
      H <- matrix(0,nrow=length(par),ncol=length(par)) 
  
     hessian.S <- function(S,ef) {
          eta <- as.numeric(mu+S[ID.coords])
            h <- units.m*exp(eta)/(1+exp(eta))
            h1 <- h/(1+exp(eta))
            q.f <- t(S)%*%R.inv%*%S
      
            grad.beta <-  t(D)%*%(y-h)

            grad.log.sigma2 <- (-n.x/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2

            grad.log.phi <- (t1.phi+0.5*as.numeric(t(S)%*%m2.phi%*%(S))/sigma2)*phi

            if(length(fixed.rel.nugget)==0) {
            	grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(S)%*%m2.nu2%*%(S))/sigma2)*nu2
            g <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2)
            } else {
            	g <- c(grad.beta,grad.log.sigma2,grad.log.phi)
            }
            
            grad2.log.lsigma2.lsigma2 <- (n.x/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
                                                    grad.log.sigma2

            grad2.log.lphi.lphi <-(t2.phi-0.5*t(S)%*%n2.phi%*%(S)/sigma2)*phi^2+
                                   grad.log.phi 
             
            if(length(fixed.rel.nugget)==0) {                       
            grad2.log.lnu2.lnu2 <- (t2.nu2-0.5*t(S)%*%n2.nu2%*%(S)/sigma2)*nu2^2+
                                         grad.log.nu2
            grad2.log.lnu2.lphi <- (t2.nu2.phi-0.5*t(S)%*%n2.nu2.phi%*%(S)/sigma2)*phi*nu2          
     
            H[1:p,1:p] <- -t(D)%*%(D*h1)
   	        H[p+1,p+1] <-  grad2.log.lsigma2.lsigma2
   	        H[p+1,p+3] <- H[p+3,p+1] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
   	        H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
   	        H[p+3,p+3] <- grad2.log.lnu2.lnu2
   	        H[p+2,p+3] <- H[p+3,p+2] <- grad2.log.lnu2.lphi
   	        H[p+2,p+2] <- grad2.log.lphi.lphi
   	        out <- list()
   	        out$mat1<- ef*(g%*%t(g)+H)                             
   	        out$g <- g*ef
   	        } else {
   	        H[1:p,1:p] <- -t(D)%*%(D*h1)
   	        H[p+1,p+1] <-  grad2.log.lsigma2.lsigma2
   	        H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
   	        H[p+2,p+2] <- grad2.log.lphi.lphi
   	        out <- list()
   	        out$mat1<- ef*(g%*%t(g)+H)                             
   	        out$g <- g*ef
   	        }
   	        out
         }            
             
         a <- rep(0,length(par))
         A <- matrix(0,length(par),length(par))
         for(i in 1:(dim(S.sim)[1])) {
   	         out.i <- hessian.S(S.sim[i,],exp.fact[i])
   	         a <- a+out.i$g
             A <- A+out.i$mat1
          }  
          (A-a%*%t(a))
   
      }
      if(messages) cat("Estimation: \n")  
      start.par <- c(par0[1:(p+1)],start.cov.pars)
      start.par[-(1:p)] <- log(start.par[-(1:p)]) 
      estim <- list()         
      if(method=="BFGS") {      	
         estimBFGS <- maxBFGS(MC.log.lik,grad.MC.log.lik,hess.MC.log.lik,
                           start.par,print.level=1*messages)  
         estim$estimate <- estimBFGS$estimate 
         estim$covariance <- solve(-estimBFGS$hessian)
         estim$log.lik <- estimBFGS$maximum                                                   
      }
      
      if(method=="nlminb") {
      	estimNLMINB <- nlminb(start.par,function(x) -MC.log.lik(x),
      	            function(x) -grad.MC.log.lik(x),
      	            function(x) -hess.MC.log.lik(x),control=list(trace=1*messages)) 
      	estim$estimate <- estimNLMINB$par                  	          
      	estim$covariance <- solve(-hess.MC.log.lik(estimNLMINB$par))
      	estim$log.lik <- -estimNLMINB$objective                 
      }
      }
      names(estim$estimate)[1:p] <- colnames(D)
      names(estim$estimate)[(p+1):(p+2)] <- c("log(sigma^2)","log(phi)")
      if(length(fixed.rel.nugget)==0) names(estim$estimate)[p+3] <- "log(nu^2)"
      rownames(estim$covariance) <- colnames(estim$covariance) <- names(estim$estimate)      
      estim$y <- y
      estim$units.m <- units.m
      estim$D <- D
      estim$coords <- coords
      estim$method <- method 
      estim$ID.coords <- ID.coords
      estim$kappa <- kappa  
      estim$h <- S.sim.res$h 
      estim$samples <- S.sim   
      if(length(fixed.rel.nugget)>0) estim$fixed.rel.nugget <- fixed.rel.nugget
      class(estim) <- "PrevMap"      
      return(estim)
}


##' @title Monte Carlo Maximum Likelihood estimation for the binomial logistic model
##' @description This function performs Monte Carlo maximum likelihood (MCML) estimation for the geostatistical binomial logistic model.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param units.m an object of class \code{\link{formula}} indicating the binomial denominators.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param ID.coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID.coords}}. These must be provided if, for example, spatial random effects are defined at household level but some of the covariates are at individual level. \bold{Warning}: the household coordinates must all be distinct otherwise see \code{\link{jitterDupCoords}}. Default is \code{NULL}.
##' @param par0 parameters of the importance sampling distribution: these should be given in the following order \code{c(beta,sigma2,phi,tau2)}, where \code{beta} are the regression coefficients, \code{sigma2} is the variance of the Gaussian process, \code{phi} is the scale parameter of the spatial correlation and \code{tau2} is the variance of the nugget effect (if included in the model).
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param kappa fixed value for the shape parameter of the Matern covariance function.
##' @param fixed.rel.nugget fixed value for the relative variance of the nugget effect; \code{fixed.rel.nugget=NULL} if this should be included in the estimation. Default is \code{fixed.rel.nugget=NULL}.
##' @param start.cov.pars a vector of length two with elements corresponding to the starting values of \code{phi} and the relative variance of the nugget effect \code{nu2}, respectively, that are used in the optimization algorithm. If \code{nu2} is fixed through \code{fixed.rel.nugget}, then \code{start.cov.pars} represents the starting value for \code{phi} only.
##' @param method method of optimization. If \code{method="BFGS"} then the \code{\link{maxBFGS}} function is used; otherwise \code{method="nlminb"} to use the \code{\link{nlminb}} function. Default is \code{method="BFGS"}.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation of the Gaussian spatial process is used when fitting the model. Default is \code{low.rank=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots that are used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the samples of the random effect is displayed after completion of conditional simulation. Default is \code{plot.correlogram=TRUE}.
##' @details
##' This function performs parameter estimation for a geostatistical binomial logistic model. Conditionally on a zero-mean stationary Gaussian process \eqn{S(x)} and mutually independent zero-mean Gaussian variables \eqn{Z} with variance \code{tau2}, the observations \code{y} are generated from a binomial distribution with probability \eqn{p} and binomial denominators \code{units.m}. A canonical logistic link is used, thus the linear predictor assumes the form
##' \deqn{\log(p/(1-p)) = d'\beta + S(x) + Z,}
##' where \eqn{d} is a vector of covariates with associated regression coefficients \eqn{\beta}. The Gaussian process \eqn{S(x)} has isotropic Matern covariance function (see \code{\link{matern}}) with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. 
##' In the \code{binomial.logistic.MCML} function, the shape parameter is treated as fixed. The relative variance of the nugget effect, \code{nu2=tau2/sigma2}, can also be fixed through the argument \code{fixed.rel.nugget}; if \code{fixed.rel.nugget=NULL}, then the relative variance of the nugget effect is also included in the estimation.
##' 
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)+Z} given the data \code{y}, in order to approximate the high-dimensiional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic epression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{maxBFGS}} (\code{method="BFGS"}), from the \pkg{maxLik} package, and \code{\link{nlminb}} (\code{method="nlminb"}).
##' 
##' \bold{Using a two-level model to include household-level and individual-level information.}
##' When analysing data from household sruveys, some of the avilable information information might be at household-level (e.g. material of house, temperature) and some at individual-level (e.g. age, gender). In this case, the Gaussian spatial process \eqn{S(x)} and the nugget effect \eqn{Z} are defined at hosuehold-level in order to account for extra-binomial variation between and within households, respectively. 
##'
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), the parameter \code{sigma2} is then multiplied by a factor \code{constant.sigma2} so as to obtain a value that is closer to the actual variance of \eqn{S(x)}. 
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: estimates of the model parameters; use the function \code{\link{coef.PrevMap}} to obtain estimates of covariance parameters on the original scale.
##' @return \code{covariance}: covariance matrix of the MCML estimates.
##' @return \code{log.lik}: maximum value of the log-likelihood.
##' @return \code{y}: binomial observations.
##' @return \code{units.m}: binomial denominators.
##' @return \code{D}: matrix of covariates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{method}: method of optimization used.
##' @return \code{ID.coords}: set of ID values defined through the argument \code{ID.coords}.
##' @return \code{kappa}: fixed value of the shape parameter of the Matern function.
##' @return \code{knots}: matrix of the spatial knots used in the low-rank approximation.
##' @return \code{const.sigma2}: adjustment factor for \code{sigma2} in the low-rank approximation.
##' @return \code{h}: vector of the values of the tuning parameter at each iteration of the Langevin-Hastings MCMC algorithm; see \code{\link{Laplace.sampling}}, or \code{\link{Laplace.sampling.lr}} if a low-rank approximation is used.
##' @return \code{samples}: matrix of the random effects samples from the importance sampling distribution used to approximate the likelihood function.
##' @return \code{fixed.rel.nugget}: fixed value for the relative variance of the nugget effect.
##' @return \code{call}: the matched call.
##' @seealso \code{\link{Laplace.sampling}}, \code{\link{Laplace.sampling.lr}}, \code{\link{summary.PrevMap}}, \code{\link{coef.PrevMap}}, \code{\link{matern}}, \code{\link{matern.kernel}},  \code{\link{control.mcmc.MCML}}, \code{\link{create.ID.coords}}.
##' @references Christensen, O. F. (2004). \emph{Monte carlo maximum likelihood in model-based geostatistics.} Journal of Computational and Graphical Statistics 13, 702-718.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @examples
##' set.seed(1234)
##' data(data_sim)
##' # Select a subset of data_sim with 50 observations
##' n.subset <- 10
##' data_subset <- data_sim[sample(1:nrow(data_sim),n.subset),]
##'
##' # Set the MCMC control parameters
##' control.mcmc <- control.mcmc.MCML(n.sim=1000,burnin=0,thin=1,
##'                            h=1.65/(n.subset^2/3))   
##'
##' # Set the parameters of the importance sampling distribution
##' par0 <- c(0,1,0.15)
##'
##' # Estimate the model parameters using MCML
##' fit.MCML <- binomial.logistic.MCML(y~1, units.m=~units.m, coords=~x1+x2,
##' data=data_subset, control.mcmc=control.mcmc,
##' fixed.rel.nugget=0,par0=par0,start.cov.pars=0.15,method="nlminb",kappa=2)            
##' summary(fit.MCML,log.cov.pars=FALSE)
##' coef(fit.MCML)
##'
##' @export
binomial.logistic.MCML <- function(formula,units.m,coords,data,ID.coords=NULL,
                                                   par0,control.mcmc,kappa,
                                                   fixed.rel.nugget=NULL,
                                                   start.cov.pars,
                                                   method="BFGS",low.rank=FALSE,
                                                   knots=NULL,
                                                   messages=TRUE,
                                                   plot.correlogram=TRUE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.") 
    if(low.rank & length(ID.coords) > 0) stop("the low-rank approximation is not available for a two-levels model with logistic link function; see instead ?binary.probit.Bayes")
    if(class(control.mcmc)!="mcmc.MCML.PrevMap") stop("control.mcmc must be of class 'mcmc.MCML.PrevMap'")
    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(class(units.m)!="formula") stop("units.m must be a 'formula' object indicating the binomial denominators.")
    if(kappa < 0) stop("kappa must be positive.")
    if(method != "BFGS" & method != "nlminb") stop("'method' must be either 'BFGS' or 'nlminb'.")
	if(!low.rank) {
		res <- geo.MCML(formula=formula,units.m=units.m,coords=coords,
		           data=data,ID.coords=ID.coords,par0=par0,control.mcmc=control.mcmc,
		           kappa=kappa,fixed.rel.nugget=fixed.rel.nugget,
		           start.cov.pars=start.cov.pars,method=method,messages=messages,
		           plot.correlogram=plot.correlogram,poisson.llik=FALSE)
	} else {
		res <- geo.MCML.lr(formula=formula,units.m=units.m,coords=coords,
		           data=data,knots=knots,par0=par0,control.mcmc=control.mcmc,
		           kappa=kappa,start.cov.pars=start.cov.pars,method=method,
		           messages=messages,plot.correlogram=plot.correlogram,poisson.llik=FALSE)		
	}
	res$call <- match.call()
	return(res)
}

##' @title Maximum Likelihood estimation for the geostatistical linear Gaussian model
##' @description This function performs maximum likelihood estimation for the geostatistical linear Gaussian Model.
##' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param ID.coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID.coords}}. These must be provided in order to define a geostatistical model where locations have multiple observations. Default is \code{ID.coords=NULL}. See the \bold{Details} section for more information.
##' @param kappa shape parameter of the Matern covariance function.
##' @param fixed.rel.nugget fixed value for the relative variance of the nugget effect; default is \code{fixed.rel.nugget=NULL} if this should be included in the estimation.
##' @param start.cov.pars if \code{ID.coords=NULL}, a vector of length two with elements corresponding to the starting values of \code{phi} and the relative variance of the nugget effect \code{nu2}, respectively, that are used in the optimization algorithm; if \code{ID.coords} is provided, a third starting value for the relative variance of the individual unexplained variation \code{nu2.star = omega2/sigma2} must be provided. If \code{nu2} is fixed through \code{fixed.rel.nugget}, then start.cov.pars represents the starting value for \code{phi} only, if \code{ID.coords=NULL}, or for \code{phi} and \code{nu2.star}, otherwise.
##' @param method method of optimization. If \code{method="BFGS"} then the \code{\link{maxBFGS}} function is used; otherwise \code{method="nlminb"} to use the \code{\link{nlminb}} function. Default is \code{method="BFGS"}.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation of the Gaussian spatial process is used when fitting the model. Default is \code{low.rank=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots that are used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param profile.llik logical; if \code{profile.llik=TRUE} the maximization of the profile likelihood is carried out. If \code{profile.llik=FALSE} the full-likelihood is used. Default is \code{profile.llik=FALSE}.
##' @details 
##' This function estimates the parameters of a geostatistical linear Gaussian model, specified as
##' \deqn{Y = d'\beta + S(x) + Z,}
##' where \eqn{Y} is the measured outcome, \eqn{d} is a vector of coavariates, \eqn{\beta} is a vector of regression coefficients, \eqn{S(x)} is a stationary Gaussian spatial process and \eqn{Z} are independent zero-mean Gaussian variables with variance \code{tau2}. More specifically, \eqn{S(x)} has an isotropic Matern covariance function with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. In the estimation, the shape parameter \code{kappa} is treated as fixed. The relative variance of the nugget effect, \code{nu2=tau2/sigma2}, can be fixed though the argument \code{fixed.rel.nugget}; if \code{fixed.rel.nugget=NULL}, then the variance of the nugget effect is also included in the estimation.
##'
##' \bold{Locations with multiple observations.}
##' If multiple observations are available at any of the sampled locations the above model is modified as follows. Let \eqn{Y_{ij}} denote the random variable associated to the measured outcome for the j-th individual at location \eqn{x_{i}}. The linear geostatistical model assumes the form \deqn{Y_{ij} = d_{ij}'\beta + S(x_{i}) + Z{i} + U_{ij},} where \eqn{S(x_{i})} and \eqn{Z_{i}} are specified as mentioned above, and \eqn{U_{ij}} are i.i.d. zer0-mean Gaussian variable with variance \eqn{\omega^2}. his model can be fitted by specifing a vector of ID for the unique set locations thourgh the argument \code{ID.coords} (see also \code{\link{create.ID.coords}}).
##'
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} can be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process, the parameter \code{sigma2} is adjusted by a factor\code{constant.sigma2}. See \code{\link{adjust.sigma2}} for more details on the the computation of the adjustment factor \code{constant.sigma2} in the low-rank approximation.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: estimates of the model parameters; use the function \code{\link{coef.PrevMap}} to obtain estimates of covariance parameters on the original scale.
##' @return \code{covariance}: covariance matrix of the ML estimates.
##' @return \code{log.lik}: maximum value of the log-likelihood.
##' @return \code{y}: response variable.
##' @return \code{D}: matrix of covariates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{ID.coords}: set of ID values defined through the argument \code{ID.coords}.
##' @return \code{method}: method of optimization used.
##' @return \code{kappa}: fixed value of the shape parameter of the Matern function.
##' @return \code{knots}: matrix of the spatial knots used in the low-rank approximation.
##' @return \code{const.sigma2}: adjustment factor for \code{sigma2} in the low-rank approximation.
##' @return \code{fixed.rel.nugget}: fixed value for the relative variance of the nugget effect.
##' @return \code{call}: the matched call.
##' @seealso \code{\link{shape.matern}}, \code{\link{summary.PrevMap}}, \code{\link{coef.PrevMap}}, \code{\link{matern}}, \code{\link{matern.kernel}}, \code{\link{maxBFGS}}, \code{\link{nlminb}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' @examples
##' data(loaloa)
##' # Empirical logit transformation
##' loaloa$logit <- log((loaloa$NO_INF+0.5)/(loaloa$NO_EXAM-loaloa$NO_INF+0.5))
##' fit.MLE <- linear.model.MLE(logit ~ 1,coords=~LONGITUDE+LATITUDE,
##'                 data=loaloa, start.cov.pars=c(0.2,0.15), 
##'                  kappa=0.5)
##' summary(fit.MLE)   
##'
##' # Low-rank approximation
##' data(data_sim)
##' n.subset <- 200
##' data_subset <- data_sim[sample(1:nrow(data_sim),n.subset),]
##'
##' # Logit transformation
##' data_subset$logit <- log(data_subset$y+0.5)/
##'                      (data_subset$units.m-
##'                       data_subset$y+0.5)                              
##' knots <- as.matrix(expand.grid(seq(-0.2,1.2,length=8),seq(-0.2,1.2,length=8)))
##' 
##' fit <- linear.model.MLE(formula=logit~1,coords=~x1+x2,data=data_subset,
##'                              kappa=2,start.cov.pars=c(0.15,0.1),low.rank=TRUE,
##'                              knots=knots)
##' summary(fit,log.cov.pars=FALSE)             
##'
linear.model.MLE <- function(formula,coords,data,ID.coords=NULL,
                             kappa,fixed.rel.nugget=NULL,
                             start.cov.pars,method="BFGS",low.rank=FALSE,
                             knots=NULL,messages=TRUE,profile.llik=FALSE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.")                                 
    if(low.rank & length(fixed.rel.nugget)>0) stop("the relative variance of the nugget effect cannot be fixed in the low-rank approximation.")     
    if(low.rank & length(ID.coords) > 0) stop("the low-rank approximation is not available for a two-levels model.")
    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(kappa < 0) stop("kappa must be positive.")
    if(length(ID.coords)==0 & profile.llik) warning("maximization of the profile likelihood is
    currently implemented only for the model with multiple observations; 
    the full-likelihood will then be used.")
    if(method != "BFGS" & method != "nlminb") stop("'method' must be either 'BFGS' or 'nlminb'.")
	if(!low.rank) {
		res <-  geo.linear.MLE(formula=formula,coords=coords,ID.coords=ID.coords,
		           data=data,kappa=kappa,fixed.rel.nugget=fixed.rel.nugget,
		           start.cov.pars=start.cov.pars,method=method,messages=messages,
                       profile.llik)
	} else {
		res <-  geo.linear.MLE.lr(formula=formula,coords=coords,
		           knots=knots,data=data,kappa=kappa,
		           start.cov.pars=start.cov.pars,method=method,messages=messages)		
	}
	res$call <- match.call()
	return(res)
}

##' @title Bayesian estimation for the geostatistical linear Gaussian model
##' @description This function performs Bayesian estimation for the geostatistical linear Gaussian model.
##' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param control.prior output from \code{\link{control.prior}}.
##' @param control.mcmc output from \code{\link{control.mcmc.Bayes}}.
##' @param kappa shape parameter of the Matern covariance function. 
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation is fitted.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @details 
##' This function performs Bayesian estimation for the geostatistical linear Gaussian model, specified as
##' \deqn{Y = d'\beta + S(x) + Z,}
##' where \eqn{Y} is the measured outcome, \eqn{d} is a vector of coavariates, \eqn{\beta} is a vector of regression coefficients, \eqn{S(x)} is a stationary Gaussian spatial process and \eqn{Z} are independent zero-mean Gaussian variables with variance \code{tau2}. More specifically, \eqn{S(x)} has an isotropic Matern covariance function with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. The shape parameter \code{kappa} is treated as fixed. 
##'
##' \bold{Priors definition.} Priors can be defined through the function \code{\link{control.prior}}. The hierarchical structure of the priors is the following. Let \eqn{\theta} be the vector of the covariance parameters \eqn{(\sigma^2,\phi,\tau^2)}; then each component of \eqn{\theta} can have independent priors freely defined by the user. However, uniform and log-normal priors are also available as default priors for each of the covariance parameters. To remove the nugget effect \eqn{Z}, no prior should be defined for \code{tau2}. Conditionally on \code{sigma2}, the vector of regression coefficients \code{beta} has a multivariate Gaussian prior with mean \code{beta.mean} and covariance matrix \code{sigma2*beta.covar}, while in the low-rank approximation the covariance matrix is simply \code{beta.covar}. 
##' 
##' \bold{Updating the covariance parameters using a Metropolis-Hastings algorithm.} In the MCMC algorithm implemented in \code{linear.model.Bayes}, the transformed parameters \deqn{(\theta_{1}, \theta_{2}, \theta_{3})=(\log(\sigma^2)/2,\log(\sigma^2/\phi^{2 \kappa}), \log(\tau^2))} are independently updated using a Metropolis Hastings algorithm. At the \eqn{i}-th iteration, a new value is proposed for each from a univariate Gaussian distrubion with variance, say \eqn{h_{i}^2}, tuned according the following adaptive scheme \deqn{h_{i} = h_{i-1}+c_{1}i^{-c_{2}}(\alpha_{i}-0.45),} where \eqn{\alpha_{i}} is the acceptance rate at the \eqn{i}-th iteration (0.45 is the optimal acceptance rate for a univariate Gaussian distribution) whilst \eqn{c_{1} > 0} and \eqn{0 < c_{2} < 1} are pre-defined constants. The starting values \eqn{h_{1}} for each of the parameters \eqn{\theta_{1}}, \eqn{\theta_{2}} and \eqn{\theta_{3}} can be set using the function \code{\link{control.mcmc.Bayes}} through the arguments \code{h.theta1}, \code{h.theta2} and \code{h.theta3}. To define values for \eqn{c_{1}} and \eqn{c_{2}}, see the documentation of \code{\link{control.mcmc.Bayes}}.
##'
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), \code{sigma2} may take very different values from the actual variance of the Gaussian process to approximate. The function \code{\link{adjust.sigma2}} can then be used to (approximately) explore the range for \code{sigma2}. For example if the variance of the Gaussian process is \code{0.5}, then an approximate value for \code{sigma2} is \code{0.5/const.sigma2}, where \code{const.sigma2} is the value obtained with \code{\link{adjust.sigma2}}.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @return An object of class "Bayes.PrevMap".
##' The function \code{\link{summary.Bayes.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: matrix of the posterior samples for each of the model parameters.
##' @return \code{S}: matrix of the posterior samplesfor each component of the random effect. This is only returned for the low-rank approximation. 
##' @return \code{y}: response variable.
##' @return \code{D}: matrix of covariarates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{kappa}: vaues of the shape parameter of the Matern function.
##' @return \code{knots}: matrix of spatial knots used in the low-rank approximation.
##' @return \code{const.sigma2}: vector of the values of the multiplicative factor used to adjust the \code{sigma2} in the low-rank approximation.
##' @return \code{h1}: vector of values taken by the tuning parameter \code{h.theta1} at each iteration.
##' @return \code{h2}: vector of values taken by the tuning parameter \code{h.theta2} at each iteration.
##' @return \code{h3}: vector of values taken by the tuning parameter \code{h.theta3} at each iteration.
##' @return \code{call}: the matched call.
##' @seealso \code{\link{control.prior}}, \code{\link{control.mcmc.Bayes}}, \code{\link{shape.matern}}, \code{\link{summary.Bayes.PrevMap}}, \code{\link{autocor.plot}}, \code{\link{trace.plot}}, \code{\link{dens.plot}}, \code{\link{matern}}, \code{\link{matern.kernel}}, \code{\link{adjust.sigma2}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' @examples
##' set.seed(1234)
##' data(loaloa)
##' # Empirical logit transformation
##' loaloa$logit <- log((loaloa$NO_INF+0.5)/(loaloa$NO_EXAM-loaloa$NO_INF+0.5))
##' 
##' cp <- control.prior(beta.mean=-2.3,beta.covar=20,
##'                              log.normal.sigma2=c(0.9,5),
##'                              log.normal.phi=c(-0.17,2),
##'                              log.normal.nugget=c(-1,1))
##' control.mcmc <- control.mcmc.Bayes(n.sim=10,burnin=0,thin=1,
##'                            h.theta1=0.5,h.theta2=0.5,h.theta3=0.5,
##'                            c1.h.theta3=0.01,c2.h.theta3=0.0001,linear.model=TRUE,
##'                            start.beta=-2.3,start.sigma2=2.45,
##'                            start.phi=0.65,start.nugget=0.34)                             
##' fit.Bayes <- linear.model.Bayes(logit ~ 1,coords=~LONGITUDE+LATITUDE,
##'                 data=loaloa,kappa=0.5, control.mcmc=control.mcmc,
##'                 control.prior = cp)
##' summary(fit.Bayes)                
linear.model.Bayes <- function(formula,coords,data,
                                                 kappa,
                                                 control.mcmc,
                                                 control.prior,
                                                 low.rank=FALSE,
                                                 knots=NULL,messages=TRUE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.")                                 
    if(low.rank & length(control.mcmc$start.nugget)==0) stop("the nugget effect must be included in the low-rank approximation.")   
    if(class(control.mcmc)!="mcmc.Bayes.PrevMap") stop("control.mcmc must be of class 'mcmc.Bayes.PrevMap'") 
    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(kappa < 0) stop("kappa must be positive.")
	if(!low.rank) {
		res <-  geo.linear.Bayes(formula=formula,coords=coords,
		           data=data,kappa=kappa,control.prior=control.prior,
		           control.mcmc=control.mcmc,messages=messages)
	} else {
		res <-  geo.linear.Bayes.lr(formula=formula,coords=coords,
		           data=data,knots=knots,kappa=kappa,
		           control.prior=control.prior,
		           control.mcmc=control.mcmc,messages=messages)		
	}
	res$call <- match.call()
	return(res)
}

##' @title Summarizing likelihood-based model fits
##' @description \code{summary} method for the class "PrevMap" that computes the standard errors and p-values of likelihood-based model fits.
##' @param object an object of class "PrevMap" obatained as result of a call to \code{\link{binomial.logistic.MCML}} or \code{\link{linear.model.MLE}}.
##' @param log.cov.pars logical; if \code{log.cov.pars=TRUE} the estimates of the covariance parameters are given on the log-scale. Note that standard errors are also adjusted accordingly. Default is \code{log.cov.pars=TRUE}.
##' @param ... further arguments passed to or from other methods.
##' @return A list with the following components
##' @return \code{linear}: logical value; \code{linear=TRUE} if a linear model was fitted and \code{linear=FALSE} otherwise.
##' @return \code{poisson}: logical value; \code{poisson=TRUE} if a Poisson model was fitted and \code{poisson=FALSE} otherwise.
##' @return \code{ck}: logical value; \code{ck=TRUE} if a low-rank approximation was used and \code{ck=FALSE} otherwise.
##' @return \code{coefficients}: matrix of the estimates, standard errors and p-values of the estimates of the regression coefficients.
##' @return \code{cov.pars}: matrix of the estimates and standard errors of the covariance parameters.
##' @return \code{log.lik}: value of likelihood function at the maximum likelihood estimates.
##' @return \code{kappa}: fixed value of the shape paramter of the Matern covariance function.
##' @return \code{fixed.rel.nugget}: fixed value for the relative variance of the nugget effect.
##' @return \code{call}: matched call.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method summary PrevMap
##' @export
summary.PrevMap <- function(object, log.cov.pars = TRUE,...) {
   res <- list()
   if(length(object$units.m)==0) {
      res$linear<-TRUE
   } else {
    	res$linear<-FALSE
   }

   if(substr(object$call[1],1,7)=="poisson") {
      res$poisson <- TRUE
   } else {
      res$poisson <- FALSE
   }
   
   if(length(dim(object$knots))==0) {
      res$ck<-FALSE
   } else {
	res$ck<-TRUE
   }	
	
   if(length(object$units.m)>0 & res$ck) {
	object$fixed.rel.nugget<-0
   }
   p <- ncol(object$D)
   object$hessian <- solve(-object$covariance)

   if(res$linear & length(object$ID.coords)>0) {
      if(length(object$fixed.rel.nugget)==0) {
         J <- diag(1,p+4)
         if(!log.cov.pars) {
            J[p+1,p+1] <- exp(object$estimate[p+1])   
            J[p+2,p+2] <- exp(object$estimate[p+2])  
            J[p+3,p+3] <- exp(object$estimate[p+3]+object$estimate[p+1])
            J[p+3,p+1] <- J[p+3,p+3]
            J[p+4,p+1] <- exp(object$estimate[p+4]+object$estimate[p+1])
            J[p+4,p+4] <- J[p+4,p+1]
            object$estimate[(p+1):(p+2)] <- exp(object$estimate[(p+1):(p+2)])
            object$estimate[(p+3):(p+4)] <- exp(object$estimate[(p+3):(p+4)])*
                                            object$estimate[p+1]
            names(object$estimate)[-(1:p)] <- c("sigma^2","phi","tau^2","omega^2")
            object$hessian <- t(J)%*%object$hessian%*%J
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- 
                                    names(object$estimate)

         } else {              
            J[p+3,p+1] <- 1
            J[p+4,p+1] <- 1
            object$estimate[(p+3):(p+4)] <- object$estimate[(p+3):(p+4)]+
                                            object$estimate[p+1]
            names(object$estimate)[-(1:p)] <- c("log(sigma^2)",
                                                "log(phi)","log(tau^2)",
                                                "log(omega^2)")
            object$hessian <- t(J)%*%object$hessian%*%J
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- 
                                    names(object$estimate)

         }
      } else {
         J <- diag(1,p+3)
         if(!log.cov.pars) {
           J[p+1,p+1] <- exp(object$estimate[p+1])   
           J[p+2,p+2] <- exp(object$estimate[p+2])  
           J[p+3,p+3] <- exp(object$estimate[p+3]+object$estimate[p+1])
           J[p+3,p+1] <- J[p+3,p+3]
           object$estimate[(p+1):(p+2)] <- exp(object$estimate[(p+1):(p+2)])
           object$estimate[p+3] <- exp(object$estimate[p+3])*
                                     object$estimate[p+1]
           names(object$estimate)[-(1:p)] <- c("sigma^2","phi","omega^2")
           object$hessian <- t(J)%*%object$hessian%*%J
           object$covar <- solve(-object$hessian)
           rownames(object$covar) <- colnames(object$covar) <- 
                                    names(object$estimate)

         } else {
           J[p+3,p+1] <- 1
           object$estimate[p+3] <- object$estimate[p+3]+
                                   object$estimate[p+1]
           names(object$estimate)[-(1:p)] <- c("log(sigma^2)",
                                               "log(phi)",
                                               "log(omega^2)")
           object$hessian <- t(J)%*%object$hessian%*%J
           object$covar <- solve(-object$hessian)
           rownames(object$covar) <- colnames(object$covar) <- 
                                    names(object$estimate)

         }
      }
   } else {
      if(length(object$fixed.rel.nugget)==0) {
         if (!log.cov.pars) {
            if(res$ck){
               res$const.sigma2 <- object$const.sigma2
        	   J <- diag(c(exp(object$estimate[p+1]),exp(object$estimate[p+2]), 
                        exp(object$estimate[p + 1] + object$estimate[p+3])))
               J[3, 1] <- exp(object$estimate[p + 1] + object$estimate[p+3])
               object$estimate[p + 1] <- exp(object$estimate[p+1])
               object$estimate[p + 2] <- exp(object$estimate[p+2])
               object$estimate[p + 3] <- object$estimate[p + 1]*exp(object$estimate[p + 3])
            } else {
        	   J <- diag(c(exp(object$estimate[(p + 1):(p + 2)]), 
                    exp(object$estimate[p + 1] + object$estimate[p+3])))
               J[3, 1] <- exp(object$estimate[p + 1] + object$estimate[p+3])
               object$estimate[p + 1] <- exp(object$estimate[p+1])
               object$estimate[p + 2] <- exp(object$estimate[p+2])
               object$estimate[p + 3] <- object$estimate[p + 1]*exp(object$estimate[p + 3])
            }
            
            names(object$estimate) <- c(names(object$estimate[1:p]), 
                "sigma^2", "phi", "tau^2")
            J <- rbind(cbind(diag(1, p), matrix(0, nrow = p, 
                ncol = 3)), cbind(matrix(0, nrow = 3, ncol = p), 
                J))
            object$hessian <- t(J) %*% object$hessian %*% J
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- 
                                    names(object$estimate)
         } else {
            if(res$ck) {
               res$const.sigma2 <- object$const.sigma2
            } 
            object$estimate[p + 3] <- object$estimate[p + 1]+object$estimate[p + 3]
            names(object$estimate) <- c(names(object$estimate[1:p]), 
                                      "log(sigma^2)", "log(phi)", "log(tau^2)")
            J <- diag(1, p + 3)
            J[p + 3, p + 1] <- 1
            object$hessian <- t(J) %*% object$hessian %*% J
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- names(object$estimate)
         }    
      } else {        
         if (!log.cov.pars) {
            if(res$ck) {
               res$const.sigma2 <- object$const.sigma2
               J <- diag(c(exp(object$estimate[p+1]),
                              exp(object$estimate[p + 2])))
               object$estimate[p + 1] <- exp(object$estimate[p +1])
               object$estimate[p + 2] <- exp(object$estimate[p + 2])
            } else {
               J <- diag(c(exp(object$estimate[p+1]),exp(object$estimate[p + 2])))
               object$estimate[p + 1] <- exp(object$estimate[p +1])
               object$estimate[p + 2] <- exp(object$estimate[p + 2])
            }   
            names(object$estimate) <- c(names(object$estimate[1:p]), 
                                      "sigma^2", "phi")
            J <- rbind(cbind(diag(1, p), matrix(0, nrow = p, 
                ncol = 2)), cbind(matrix(0, nrow = 2, ncol = p), 
                J))
            object$hessian <- t(J) %*% object$hessian %*% J
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- names(object$estimate)
         } else {
            if(res$ck) {
               res$const.sigma2 <- object$const.sigma2	
            } 
            names(object$estimate) <- c(names(object$estimate[1:p]), 
                                     "log(sigma^2)", "log(phi)")
            object$covar <- solve(-object$hessian)
            rownames(object$covar) <- colnames(object$covar) <- names(object$estimate)
         }
      }
   }   
   se <- sqrt(diag(object$covar))
   zval <- object$estimate[1:p]/se[1:p]
   TAB <- cbind(Estimate = object$estimate[1:p], StdErr = se[1:p], 
   z.value = zval, p.value = 2 * pnorm(-abs(zval)))
   cov.pars <- cbind(Estimate = object$estimate[-(1:p)],StdErr = se[-(1:p)])
   res$coefficients <- TAB
   res$cov.pars <- cov.pars     
   res$log.lik <- object$log.lik
   res$kappa <- object$kappa
   res$fixed.rel.nugget <- object$fixed.rel.nugget 
   res$call <- object$call               
   class(res) <- "summary.PrevMap"
   return(res)
}

##' @title Trace-plots of the importance sampling distribution samples from the MCML method 
##' @description Trace-plots of the MCMC samples from the importance sampling distribution used in \code{\link{binomial.logistic.MCML}}.
##' @param object an object of class "PrevMap" obatained as result of a call to \code{\link{binomial.logistic.MCML}}.
##' @param component a positive integer indicating the number of the random effect component for which a trace-plot is required. If \code{component=NULL}, then a component is selected at random. Default is \code{component=NULL}. 
##' @param ... further arguments passed to \code{\link{plot}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
trace.plot.MCML <- function(object,component=NULL,...) {
	if(class(object)!="PrevMap") {
	   stop("'object' must be of class 'PrevMap'")	
	}
	
	if(length(object$units.m)==0) {
	   stop("'object' must be the output of a call to the 'binomial.logistic.MCML' function")
	}
	
	n.samples <- ncol(object$samples)
	if(length(component)==0) {
		component <- sample(1:n.samples,1)
	}
	
	if(component > n.samples) stop("'components' must not be greater than the number of random effects in the model")
	
	re <- object$samples[,component]
	cat("Plotted component number",component,"\n")
	plot(re,type="l",ylab="",xlab="Iteration",
	     main=paste("Component number",component))
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method print summary.PrevMap
##' @export
print.summary.PrevMap <- function(x,...) {
   if(x$ck==FALSE) {	
      if(x$linear) {
   	   cat("Geostatistical linear model \n")
      } else if(x$poisson) {
         cat("Geostatistical Poisson model \n")
      } else {
         cat("Geostatistical binomial model \n")
      }  
      cat("Call: \n")
      print(x$call)
      cat("\n")	
      printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
      cat("\n")
      if(x$linear) {
   	   cat("Log-likelihood: ",x$log.lik,"\n \n",sep="")
      } else {
   	   cat("Objective function: ",x$log.lik,"\n \n",sep="")
      }	 
      if(length(x$fixed.rel.nugget)==0) {
         cat("Covariance parameters Matern function (kappa=",x$kappa,") \n",sep="")
         printCoefmat(x$cov.pars,P.values=FALSE)
      } else {
         cat("Covariance parameters Matern function \n")
         cat("(fixed relative variance tau^2/sigma^2= ",x$fixed.rel.nugget,") \n",sep="")
         printCoefmat(x$cov.pars,P.values=FALSE)	
      }
   } else {
      if(x$linear) {
         cat("Geostatistical linear model \n")
      } else if(x$poisson) {
         cat("Geostatistical Poisson model \n")
      } else {
         cat("Geostatistical binomial model \n")
      } 	
      cat("(low-rank approximation) \n")
      cat("Call: \n")
      print(x$call)
      cat("\n")	
      printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
      cat("\n")
      if(x$linear) {
   	   cat("Log-likelihood: ",x$log.lik,"\n \n",sep="")
      } else {
   	   cat("Objective function: ",x$log.lik,"\n \n",sep="")
      }	 
      cat("Matern kernel parameters (kappa=",x$kappa,") \n",sep="")
      cat("Adjustment factorfor sigma^2: ",x$const.sigma2 ,"\n",sep="")
      printCoefmat(x$cov.pars,P.values=FALSE)   	
   }
   cat("\n")
   cat("Legend: \n")
   cat("sigma^2 = variance of the Gaussian process \n")
   cat("phi = scale of the spatial correlation \n")
   if(length(x$fixed.rel.nugget)==0) cat("tau^2 = variance of the nugget effect \n")
   if(any(names(x$cov.pars[,1])=="omega^2") | 
      any(names(x$cov.pars[,1])=="log(omega^2)")) {
      cat("omega^2 = variance of the individual unexplained variation \n")
   }
}

##' @title Spatial predictions for the binomial logistic model using plug-in of MCML estimates
##' @description This function performs spatial prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a geostatistical binomial logistic model.
##' @param object an object of class "PrevMap" obtained as result of a call to \code{\link{binomial.logistic.MCML}}.
##' @param grid.pred a matrix of prediction locations.
##' @param predictors a data frame of the values of the explanatory variables at each of the locations in \code{grid.pred}; each column correspond to a variable and each row to a location. \bold{Warning:} the names of the columns in the data frame must match those in the data used to fit the model. Default is \code{predictors=NULL} for models with only an intercept.
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param type a character indicating the type of spatial predictions: \code{type="marginal"} for marginal predictions or \code{type="joint"} for joint predictions. Default is \code{type="marginal"}. In the case of a low-rank approximation only joint predictions are available.
##' @param scale.predictions a character vector of maximum length 3, indicating the required scale on which spatial prediction is carried out: "logit", "prevalence" and "odds". Default is \code{scale.predictions=c("logit","prevalence","odds")}.
##' @param quantiles a vector of quantiles used to summarise the spatial predictions.
##' @param standard.errors logical; if \code{standard.errors=TRUE}, then standard errors for each \code{scale.predictions} are returned. Default is \code{standard.errors=FALSE}.
##' @param thresholds a vector of exceedance thresholds; default is \code{thresholds=NULL}.
##' @param scale.thresholds a character value indicating the scale on which exceedance thresholds are provided; \code{"logit"}, \code{"prevalence"} or \code{"odds"}. Default is \code{scale.thresholds=NULL}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the conditional simulations is displayed. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return A "pred.PrevMap" object list with the following components: \code{logit}; \code{prevalence}; \code{odds}; \code{exceedance.prob}, corresponding to a matrix of the exceedance probabilities where each column corresponds to a specified value in \code{thresholds}; \code{samples}, corresponding to a matrix of the predictive samples at each prediction locations for the linear predictor of the binomial logistic model (if \code{scale.predictions="logit"} this component is \code{NULL}); \code{grid.pred} prediction locations. 
##' Each of the three components \code{logit}, \code{prevalence} and  \code{odds} is also a list with the following components:
##' @return \code{predictions}: a vector of the predictive mean for the associated quantity (logit, odds or prevalence).
##' @return \code{standard.errors}: a vector of prediction standard errors (if \code{standard.errors=TRUE}).
##' @return \code{quantiles}: a matrix of quantiles of the resulting predictions with each column corresponding to a quantile specified through the argument \code{quantiles}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom geoR matern 
##' @export 
spatial.pred.binomial.MCML <- function(object,grid.pred,predictors=NULL,control.mcmc,
                                                    type="marginal",
                                                    scale.predictions=c("logit","prevalence","odds"),
                                                    quantiles=c(0.025,0.975),
                                                    standard.errors=FALSE,                                                                                                                                                      
                                                    thresholds=NULL,
                                                    scale.thresholds=NULL,
                                                    plot.correlogram=FALSE,
                                                    messages=TRUE) {
    if(nrow(grid.pred) < 2) stop("prediction locations must be at least two.")       
    if(length(predictors)>0 && class(predictors)!="data.frame") stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.")
    if(length(predictors)>0 && any(is.na(predictors))) stop("missing values found in 'predictors'.")
	p <- object$p <- ncol(object$D)	
	kappa <- object$kappa	
	n.pred <- nrow(grid.pred)
	coords <- object$coords
    if(type=="marginal" & length(object$knots) > 0) warning("only joint predictions are avilable for the low-rank approximation.")
    if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions should be marginal or joint")
	for(i in 1:length(scale.predictions)) {
		if(any(c("logit","prevalence","odds")==scale.predictions[i])==
		          FALSE) stop("invalid scale.predictions.")
	}
	
	if(length(thresholds)>0) {
		if(any(scale.predictions==scale.thresholds)==FALSE) {
			stop("scale thresholds must be equal to a scale prediction")
		}
	}
	
	if(length(thresholds)==0 & length(scale.thresholds)>0 |
       length(thresholds)>0 & length(scale.thresholds)==0) stop("to estimate exceedance probabilities both thresholds and scale.thresholds.")

	if(object$p==1) {
	   predictors <- matrix(1,nrow=n.pred)	
	} else {
	   if(length(dim(predictors))==0) stop("covariates at prediction locations should be provided.")	
	   predictors <- as.matrix(model.matrix(delete.response(terms(formula(object$call))),data=predictors))
	   if(nrow(predictors)!=nrow(grid.pred)) stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")
	   if(ncol(predictors)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
	}
	
	if(length(dim(object$knots)) > 0) {	
	beta <- object$estimate[1:p]
	sigma2 <- exp(object$estimate["log(sigma^2)"])/object$const.sigma2
	rho <- exp(object$estimate["log(phi)"])*2*sqrt(object$kappa)
	knots <- object$knots
	U.k <- as.matrix(pdist(coords,knots))
	K <- matern.kernel(U.k,rho,kappa)
	mu.pred <- as.numeric(predictors%*%beta)
	object$mu <- object$D%*%beta	
	Z.sim.res <- Laplace.sampling.lr(object$mu,sigma2,K,
	object$y,object$units.m,control.mcmc,
	plot.correlogram=plot.correlogram,messages=messages,poisson.llik=FALSE)
	Z.sim <- Z.sim.res$samples
	U.k.pred <- as.matrix(pdist(grid.pred,knots))
	K.pred <- matern.kernel(U.k.pred,rho,kappa) 
	out <- list()     
    eta.sim <- sapply(1:(dim(Z.sim)[1]), function(i) mu.pred+K.pred%*%Z.sim[i,])           
    
    if(any(scale.predictions=="logit")) {
    	   if(messages) cat("Spatial predictions: logit \n")    	   
    	   out$logit$predictions <- apply(eta.sim,1,mean)
       if(standard.errors) {
          out$logit$standard.errors <- apply(eta.sim,1,sd)
       }
       if(length(quantiles) > 0) {
          out$logit$quantiles <- t(apply(eta.sim,1,function(r) quantile(r,quantiles)))         
       }
       
       if(length(thresholds) > 0 && scale.thresholds=="logit") {
       	  out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	  colnames(out$exceedance.prob) <- paste(thresholds,sep="")
          for(j in 1:length(thresholds)) {
             out$exceedance.prob[,j] <- apply(eta.sim,1,function(r) mean(r > thresholds[j]))	
          } 
       }		   
    }
    if(any(scale.predictions=="odds")) {
       if(messages) cat("Spatial predictions: odds \n") 
       odds.sim <- exp(eta.sim)    	   
       out$odds$predictions <- apply(odds.sim,1,mean)
       if(standard.errors) {
          out$odds$standard.errors <- apply(odds.sim,1,sd)
       }
       
       if(length(quantiles) > 0) {
          out$odds$quantiles <- t(apply(odds.sim,1,function(r) quantile(r,quantiles)))         
       }
       	
       if(length(thresholds) > 0 && scale.thresholds=="odds") {
       	  out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	  colnames(out$exceedance.prob) <- paste(thresholds,sep="")
          for(j in 1:length(thresholds)) {
             out$exceedance.prob[,j] <- apply(odds.sim,1,function(r) mean(r > thresholds[j]))	
          } 
       }	   
    }  
	} else {	
	beta <- object$estimate[1:p]
	sigma2 <- exp(object$estimate[p+1])
	phi <- exp(object$estimate[p+2])
	if(length(object$fixed.rel.nugget)==0){		
	   tau2 <- sigma2*exp(object$estimate[p+3]) 
	} else {
	   tau2 <- object$fixed.rel.nugget*sigma2
	}
	

	U <- dist(coords)
	U.pred.coords <- as.matrix(pdist(grid.pred,coords))	
	Sigma <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov 
	Sigma.inv <- solve(Sigma)   	      
    C <- sigma2*matern(U.pred.coords,phi,kappa)	                
    A <- C%*%Sigma.inv
	out <- list()  
	mu.pred <- as.numeric(predictors%*%beta)   
	object$mu <- object$D%*%beta      	             	
    if(length(object$ID.coords)>0) {
    	   S.sim.res <- Laplace.sampling(object$mu,Sigma,
	              object$y,object$units.m,control.mcmc,
	              object$ID.coords,
	              plot.correlogram=plot.correlogram,messages=messages)
	       S.sim <- S.sim.res$samples      
       mu.cond <- sapply(1:(dim(S.sim)[1]),function(i) mu.pred+A%*%S.sim[i,]) 	              
    } else {
	   S.sim.res <- Laplace.sampling(object$mu,Sigma,
	                 object$y,object$units.m,control.mcmc,
	                 plot.correlogram=plot.correlogram,messages=messages) 	   	
       	S.sim <- S.sim.res$samples     
       mu.cond <- sapply(1:(dim(S.sim)[1]),function(i) mu.pred+A%*%(S.sim[i,]-object$mu))   	
    }
	              
    if(type=="marginal") {
    	   if(messages) cat("Type of predictions:",type,"\n")
    	   sd.cond <- sqrt(sigma2-diag(A%*%t(C)))
    } else if (type=="joint") {
    	   if(messages) cat("Type of predictions: ",type," (this step might be demanding) \n")
    	   Sigma.pred <-  varcov.spatial(coords=grid.pred,cov.model="matern",
	                cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
	   Sigma.cond <- Sigma.pred - A%*%t(C)
	   sd.cond <- sqrt(diag(Sigma.cond))	   
    }  
    
    if((length(quantiles) > 0) | 
       (any(scale.predictions=="prevalence")) | 
       (any(scale.predictions=="odds")) |
       (length(scale.thresholds)>0)) {
       if(type=="marginal") {
          eta.sim <- sapply(1:(dim(S.sim)[1]), function(i) rnorm(n.pred,mu.cond[,i],sd.cond))
       } else if(type=="joint") {
          Sigma.cond.sroot <- t(chol(Sigma.cond))
          eta.sim <- sapply(1:(dim(S.sim)[1]), function(i) mu.cond[,i]+
                                                                   Sigma.cond.sroot%*%rnorm(n.pred))                                                                   
       }      	
    }
    
        if(any(scale.predictions=="logit")) {
    	   if(messages) cat("Spatial predictions: logit \n")    	   
    	   out$logit$predictions <- apply(mu.cond,1,mean)
       if(standard.errors) {
          out$logit$standard.errors <- sqrt(sd.cond^2+diag(A%*%cov(S.sim)%*%t(A)))
       }
       if(length(quantiles) > 0) {
          out$logit$quantiles <- t(apply(eta.sim,1,function(r) quantile(r,quantiles)))         
       }
       
       if(length(thresholds) > 0 && scale.thresholds=="logit") {
       	  out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	  colnames(out$exceedance.prob) <- paste(thresholds,sep="")
          for(j in 1:length(thresholds)) {
             out$exceedance.prob[,j] <- apply(eta.sim,1,function(r) mean(r > thresholds[j]))	
          } 
       }		   
    }
    
    if(any(scale.predictions=="odds")) {
       if(messages) cat("Spatial predictions: odds \n") 
       odds.sim <- exp(eta.sim)    	   
       out$odds$predictions <- apply(exp(mu.cond+0.5*sd.cond^2),1,mean)
       if(standard.errors) {
          out$odds$standard.errors <- apply(odds.sim,1,sd)
       }
       
       if(length(quantiles) > 0) {
          out$odds$quantiles <- t(apply(odds.sim,1,function(r) quantile(r,quantiles)))         
       }
       	
       if(length(thresholds) > 0 && scale.thresholds=="odds") {
       	  out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	  colnames(out$exceedance.prob) <- paste(thresholds,sep="")
          for(j in 1:length(thresholds)) {
             out$exceedance.prob[,j] <- apply(odds.sim,1,function(r) mean(r > thresholds[j]))	
          } 
       }	   
    }  	
    }
       
    
    if(any(scale.predictions=="prevalence")) {
       if(messages) cat("Spatial predictions: prevalence \n") 
       prev.sim <- exp(eta.sim)/(1+exp(eta.sim))    	   
       out$prevalence$predictions <- apply(prev.sim,1,mean)
       if(standard.errors) {
          out$prevalence$standard.errors <- apply(prev.sim,1,sd)
       }
       
       if(length(quantiles) > 0) {
          out$prevalence$quantiles <- t(apply(prev.sim,1,function(r) quantile(r,quantiles)))         
       }
       	
       if(length(thresholds) > 0 && scale.thresholds=="prevalence") {       	
       	  out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	  colnames(out$exceedance.prob) <- paste(thresholds,sep="")
          for(j in 1:length(thresholds)) {
             out$exceedance.prob[,j] <- apply(prev.sim,1,function(r) mean(r > thresholds[j]))	
          } 
       }	   
    }
    
    if(any(scale.predictions=="odds") | 
       any(scale.predictions=="prevalence")) {
       out$samples <- eta.sim	
    }
    out$grid <- grid.pred
    class(out) <- "pred.PrevMap"
    out
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
##' @importFrom pdist pdist
geo.MCML.lr <- function(formula,units.m,coords,data,knots,
                                                   par0,control.mcmc,kappa,
                                                   start.cov.pars,
                                                   method,plot.correlogram,messages,
                                                   poisson.llik) {
    
    knots <- as.matrix(knots)
    start.cov.pars <- as.numeric(start.cov.pars)
    if(length(start.cov.pars)!=1) stop("wrong values for start.cov.pars")
    kappa <- as.numeric(kappa)
    coords <- as.matrix(model.frame(coords,data))
       
    if(any(is.na(data))) stop("missing values are not accepted")
    if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")                                                        
    if(any(method==c("BFGS","nlminb"))==FALSE) stop("method must be either BFGS or nlminb.")
   
    der.rho <- function(u,rho,kappa) {
        u <- u+10e-16
        if(kappa==2) {
           out <- ((2^(7/2)*u-4*rho)*exp(-(2^(3/2)*u)/rho))/(sqrt(pi)*rho^3)
        } else {
           out <- (kappa^(kappa/4+1/4)*sqrt(gamma(kappa+1))*rho^(-kappa/2-5/2)*u^(kappa/2-1/2)*
                     (2*sqrt(kappa)*besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*
                     u+2*besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*sqrt(kappa)*
                     u-besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*rho-
                     besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*rho))/
                     (sqrt(gamma(kappa))*gamma((kappa+1)/2)*sqrt(pi))
        }
        out
     }

    der2.rho <- function(u,rho,kappa) {
        u <- u+10e-16
        if(kappa==2) {
          	out <- (8*(4*u^2-2^(5/2)*rho*u+rho^2)*exp(-(2^(3/2)*u)/rho))/(sqrt(pi)*rho^5)
        } else {
            out <- (kappa^(kappa/4+1/4)*sqrt(gamma(kappa+1))*rho^(-kappa/2-9/2)*u^(kappa/2-1/2)*
                      (4*kappa*besselK((2*sqrt(kappa)*u)/rho,(kappa+3)/2)*u^2+
                      8*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*u^2+
                      4*besselK((2*sqrt(kappa)*u)/rho,(kappa-5)/2)*kappa*u^2-
                      4*kappa^(3/2)*besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*rho*u-12*sqrt(kappa)*
                      besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*rho*u-4*
                      besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*kappa^(3/2)*rho*u-
                      12*besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*sqrt(kappa)*rho*u+
                      besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa^2*rho^2+
                      4*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*rho^2+
                      3*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*rho^2))/
                     (2*sqrt(gamma(kappa))*gamma((kappa+1)/2)*sqrt(pi))
        }
        out
    }
    mf <- model.frame(formula,data=data)
	y <- as.numeric(model.response(mf))
	n <- length(y)
    if(poisson.llik && length(units.m)==0) {
      units.m <- rep(1,n)
    } else {
      units.m <-  as.numeric(model.frame(units.m,data)[,1])  
    }
    if(is.numeric(units.m)==FALSE) stop("'units.m' must be numeric.")
    
	D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
	p <- ncol(D)
	if(any(par0[-(1:p)] <= 0)) stop("the covariance parameters in 'par0' must be positive.")
	beta0 <- par0[1:p]
	mu0 <- as.numeric(D%*%beta0)
		
	sigma2.0 <- par0[p+1]
	rho0 <- par0[p+2]*2*sqrt(kappa)	
	
	U <- as.matrix(pdist(coords,knots))
	N <- nrow(knots)
	K0 <- matern.kernel(U,rho0,kappa)
      const.sigma2.0 <- mean(apply(K0,1,function(r) sqrt(sum(r^2))))
      sigma2.0 <- sigma2.0/const.sigma2.0
    
      n.sim <- control.mcmc$n.sim
	Z.sim.res <- Laplace.sampling.lr(mu0,sigma2.0,K0,y,units.m,control.mcmc,
	              plot.correlogram=plot.correlogram,
                  messages=messages,poisson.llik=poisson.llik)
      Z.sim <- Z.sim.res$samples
	log.integrand <- function(Z,val,K) {
         S <- as.numeric(K%*%Z)
         eta <- as.numeric(val$mu+S)

         if(poisson.llik) {
            llik <- sum(y*eta-units.m*exp(eta))
         } else {
            llik <- sum(y*eta-units.m*log(1+exp(eta)))
         }
	   
         -0.5*(N*log(val$sigma2)+sum(Z^2)/val$sigma2)+
         llik        
      }    
   
      compute.log.f <- function(sub.par,K) {
         beta <- sub.par[1:p];
         val <- list()
         val$sigma2 <- exp(sub.par[p+1])
         val$mu <- as.numeric(D%*%beta)
         sapply(1:(dim(Z.sim)[1]),function(i) log.integrand(Z.sim[i,],val,K))
      }            

      log.f.tilde <- compute.log.f(c(beta0,log(sigma2.0)),K0)

      MC.log.lik <- function(par) {
         rho <- exp(par[p+2])
         sub.par <- par[-(p+2)]	
         K <- matern.kernel(U,rho,kappa)
         log(mean(exp(compute.log.f(sub.par,K)-log.f.tilde)))
      }          
    
      grad.MC.log.lik <- function(par) {
         beta <- par[1:p]; mu <- D%*%beta
         sigma2 <- exp(par[p+1])
         rho <- exp(par[p+2]) 

         K <- matern.kernel(U,rho,kappa)                           
                             
         exp.fact <- exp(compute.log.f(par[-(p+2)],K)-log.f.tilde)
         L.m <- sum(exp.fact)
         exp.fact <- exp.fact/L.m           
        
         D.rho <- der.rho(U,rho,kappa)                        
                                    
         gradient.Z <- function(Z) {
            S <- as.numeric(K%*%Z)
            eta <- as.numeric(mu+S)
            if(poisson.llik) {
               h <- units.m*exp(eta)
            } else {
               h <- units.m*exp(eta)/(1+exp(eta))
            }
            grad.beta <- as.numeric(t(D)%*%(y-h))
            S.rho <- as.numeric(D.rho%*%Z)
           
            grad.log.sigma2 <- -0.5*(N/sigma2-sum(Z^2)/(sigma2^2))*sigma2         

            grad.log.rho <- (t(S.rho)%*%(y-h))*rho

            c(grad.beta,grad.log.sigma2,grad.log.rho)
         } 
         out <- rep(0,length(par))
         for(i in 1:(dim(Z.sim)[1])) {
   	        out <- out + exp.fact[i]*gradient.Z(Z.sim[i,])
         }
         out
      }
   
      hess.MC.log.lik <- function(par) {
        beta <- par[1:p]; mu <- D%*%beta
        sigma2 <- exp(par[p+1])
        rho <- exp(par[p+2]) 

        K <- matern.kernel(U,rho,kappa)                           
                             
        exp.fact <- exp(compute.log.f(par[-(p+2)],K)-log.f.tilde)
        L.m <- sum(exp.fact)
        exp.fact <- exp.fact/L.m           
        
        D1.rho <- der.rho(U,rho,kappa)
        D2.rho <- der2.rho(U,rho,kappa)
        H <- matrix(0,nrow=length(par),ncol=length(par)) 
        
        hessian.Z <- function(Z,ef) {
           S <- as.numeric(K%*%Z)
           eta <- as.numeric(mu+S)
           
           if(poisson.llik) {
              h <- units.m*exp(eta)
           } else {
              h <- units.m*exp(eta)/(1+exp(eta))
           }

           grad.beta <- as.numeric(t(D)%*%(y-h))
           S1.rho <- as.numeric(D1.rho%*%Z)
           S2.rho <- as.numeric(D2.rho%*%Z)          

           grad.log.sigma2 <- -0.5*(N/sigma2-sum(Z^2)/(sigma2^2))*sigma2         

           grad.log.rho <- (t(S1.rho)%*%(y-h))*rho

           g <- c(grad.beta,grad.log.sigma2,grad.log.rho) 
           
           if(poisson.llik) {
              h1 <- units.m*exp(eta)
           } else {
              h1 <- units.m*exp(eta)/((1+exp(eta))^2)
           }
           
           grad2.log.beta.beta <- -t(D)%*%(D*h1)      
           grad2.log.beta.lrho <- -t(D)%*%(as.vector(D1.rho%*%Z)*h1)*rho                      
           
     
           grad2.log.lsigma2.lsigma2 <- -0.5*(-N/(sigma2^2)+2*sum(Z^2)/(sigma2^3))*sigma2^2+
                                                          grad.log.sigma2
   
           grad2.log.lrho.lrho <- (t(S2.rho)%*%(y-h)-t((S1.rho)^2)%*%h1)*rho^2+
                                             grad.log.rho         

           H[1:p,1:p] <-  grad2.log.beta.beta
           H[1:p,p+2] <- H[p+2,1:p]<- grad2.log.beta.lrho
   	     H[p+1,p+1] <-  grad2.log.lsigma2.lsigma2
   	     H[p+2,p+2] <- grad2.log.lrho.lrho
                                                       
     
   	     out <- list()
   	     out$mat1<- ef*(g%*%t(g)+H)                             
   	     out$g <- g*ef
   	     out
        }            
             
        a <- rep(0,length(par))
        A <- matrix(0,length(par),length(par))
        for(i in 1:(dim(Z.sim)[1])) {
   	       out.i <- hessian.Z(Z.sim[i,],exp.fact[i])
   	       a <- a+out.i$g
           A <- A+out.i$mat1
        }  
        (A-a%*%t(a))
      }

      if(messages) cat("Estimation: \n")
      start.par <- rep(NA,p+2)
      start.par[1:p] <- par0[1:p]
      start.par[p+1] <- log(sigma2.0)
      start.par[p+2] <- log(2*sqrt(kappa)*start.cov.pars)
      
      estim <- list()
      if(method=="BFGS") {
         estimBFGS <- maxBFGS(MC.log.lik,grad.MC.log.lik,hess.MC.log.lik,
                           start.par,print.level=1*messages)   
         estim$estimate <- estimBFGS$estimate  
         estim$covariance <- matrix(NA,nrow=p+2,ncol=p+2)
         hessian <- estimBFGS$hessian
         estim$log.lik <- estimBFGS$maximum                                                 
      }
      
      if(method=="nlminb") {
      	estimNLMINB <- nlminb(start.par,function(x) -MC.log.lik(x),
      	            function(x) -grad.MC.log.lik(x),
      	            function(x) -hess.MC.log.lik(x),control=list(trace=1*messages)) 
      	estim$estimate <- estimNLMINB$par            
      	estim$covariance <- matrix(NA,nrow=p+2,ncol=p+2)
      	hessian <- hess.MC.log.lik(estimNLMINB$par)
      	estim$log.lik <- -estimNLMINB$objective              
      }
      K.hat <- matern.kernel(U,exp(estim$estimate[p+2]),kappa)
      const.sigma2 <- mean(apply(K.hat,1,function(r) sqrt(sum(r^2))))
      const.sigma2.grad <- mean(apply(U,1,function(r) {
      	                          rho <- exp(estim$estimate[p+2])
		                          k <- matern.kernel(r,rho,kappa)
		                          k.grad <- der.rho(r,rho,kappa)
		                          den <- sqrt(sum(k^2))
		                          num <- sum(k*k.grad)
	                              num/den
	                              }))
      estim$estimate[p+1] <- estim$estimate[p+1]+log(const.sigma2)
      estim$estimate[p+2] <- estim$estimate[p+2]-log(2*sqrt(kappa))      
      J <- diag(1,p+2)
      J[p+1,p+2] <- exp(estim$estimate[p+2])*const.sigma2.grad/const.sigma2
      hessian <- J%*%hessian%*%t(J) 
      names(estim$estimate)[1:p] <- colnames(D)
      names(estim$estimate)[(p+1):(p+2)] <- c("log(sigma^2)","log(phi)") 
      estim$covariance <- solve(-J%*%hessian%*%t(J))
      colnames(estim$covariance) <- rownames(estim$covariance) <- 
      names(estim$estimate)
      estim$y <- y
      estim$units.m <- units.m
      estim$D <- D
      estim$coords <- coords
      estim$method <- method 
      estim$kappa <- kappa
      estim$knots <- knots    
      estim$const.sigma2 <- const.sigma2      
      estim$h <- Z.sim.res$h
      estim$samples <- Z.sim
      class(estim) <- "PrevMap"
      return(estim)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom geoR varcov.spatial
##' @importFrom maxLik maxBFGS
geo.linear.MLE <- function(formula,coords,data,ID.coords,
                  kappa,fixed.rel.nugget=NULL,start.cov.pars,
                  method="BFGS",messages=TRUE,profile.llik) {              
   start.cov.pars <- as.numeric(start.cov.pars)
   if(any(start.cov.pars<0)) stop("start.cov.pars must be positive.")   
   kappa <- as.numeric(kappa)
   if(any(method==c("BFGS","nlminb"))==FALSE) stop("method must be either BFGS or nlminb.")                                       	                                     	
   mf <- model.frame(formula,data=data)
   y <- as.numeric(model.response(mf))
   if(length(ID.coords)>0) {
      m <- length(y)
   } else {
      n <- length(y)
   }
   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
   coords <- as.matrix(model.frame(coords,data))
   if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("wrong set of coordinates.")
   p <- ncol(D)
   if(length(fixed.rel.nugget)>0){
      if(length(fixed.rel.nugget) != 1 | fixed.rel.nugget < 0) stop("negative fixed nugget value or wrong length of fixed.rel.nugget")
      if(length(ID.coords)>0) {
         if(length(start.cov.pars)!=2) stop("wrong length of start.cov.pars")	
      } else {
         if(length(start.cov.pars)!=1) stop("wrong length of start.cov.pars")	
      }
      if(messages) cat("Fixed relative variance of the nugget effect:",fixed.rel.nugget,"\n")	
   } else {
	if(length(ID.coords)>0) {
         if(length(start.cov.pars)!=3) stop("wrong length of start.cov.pars")	
      } else {
         if(length(start.cov.pars)!=2) stop("wrong length of start.cov.pars")	
      }		
   }	
   
   if(length(ID.coords)>0) coords <- unique(coords)
   U <- dist(coords)
   
   if(length(ID.coords)==0) {
      if(length(fixed.rel.nugget)==0) {
         R <- varcov.spatial(dists.lowertri=U,cov.model="matern",kappa=kappa,
    	                               cov.pars=c(1,start.cov.pars[1]),nugget=start.cov.pars[2])$varcov
      } else {
         R <- varcov.spatial(dists.lowertri=U,cov.model="matern",kappa=kappa,
    	                               cov.pars=c(1,start.cov.pars[1]),nugget=fixed.rel.nugget)$varcov
      }	
                                   
      R.inv <- solve(R)
   
      beta.start <- as.numeric(solve(t(D)%*%R.inv%*%D)%*%t(D)%*%R.inv%*%y)
      diff.b <- y-D%*%beta.start
      sigma2.start <- as.numeric(t(diff.b)%*%R.inv%*%diff.b/length(y))
      start.par <- c(beta.start,sigma2.start,start.cov.pars)                           
    
      start.par[-(1:p)] <- log(start.par[-(1:p)])
   } else {
      
      n.coords <- as.numeric(tapply(ID.coords,ID.coords,length))
      DtD <- t(D)%*%D
      Dty <- as.numeric(t(D)%*%y)
      D.tilde <- sapply(1:p,function(i) as.numeric(tapply(D[,i],ID.coords,sum)))
      y.tilde <- tapply(y,ID.coords,sum)

      if(profile.llik) {
         start.par <- log(start.cov.pars)
      } else {
         
         phi.start <- start.cov.pars[1]                  
         if(length(fixed.rel.nugget)==1) {
            nu2.1.start <- fixed.rel.nugget
            nu2.2.start <- start.cov.pars[2]
         } else {
            nu2.1.start <- start.cov.pars[2]
            nu2.2.start <- start.cov.pars[3]
         }
         R.start <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi.start),
                 kappa=kappa)$varcov
         diag(R.start) <- diag(R.start)+nu2.1.start
         R.start.inv <- solve(R.start)
         Omega <- R.start.inv 
         diag(Omega) <- (diag(Omega)+n.coords/nu2.2.start)

         Omega.inv <- solve(Omega)
         M.beta <- DtD/nu2.2.start+
                -t(D.tilde)%*%Omega.inv%*%D.tilde/(nu2.2.start^2)
         v.beta <- Dty/nu2.2.start-t(D.tilde)%*%Omega.inv%*%y.tilde/
                   (nu2.2.start^2)
         beta.start <- as.numeric(solve(M.beta)%*%v.beta)
      
         M.det <- t(R.start*n.coords/nu2.2.start)
         diag(M.det) <- diag(M.det)+1
      
         mu.start <- as.numeric(D%*%beta.start)
         diff.y <- y-mu.start 
         diff.y.tilde <- tapply(diff.y,ID.coords,sum)
         sigma2.start <- as.numeric((sum(diff.y^2)/nu2.2.start-
                       t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde/
                       (nu2.2.start^2))/m)
     
         if(length(fixed.rel.nugget)==0) {
            start.par <- c(beta.start,log(c(sigma2.start,phi.start,
                                            nu2.1.start,nu2.2.start)))
         } else {
            start.par <- c(beta.start,log(c(sigma2.start,phi.start,
                                            nu2.2.start)))
         }
      }
   } 

   der.phi <- function(u,phi,kappa) {
      u <- u+10e-16
      if(kappa==0.5) {
         out <- (u*exp(-u/phi))/phi^2
      } else {
        out <- ((besselK(u/phi,kappa+1)+besselK(u/phi,kappa-1))*
                  phi^(-kappa-2)*u^(kappa+1))/(2^kappa*gamma(kappa))-
                  (kappa*2^(1-kappa)*besselK(u/phi,kappa)*phi^(-kappa-1)*
                  u^kappa)/gamma(kappa)
      }
      out
   }

   der2.phi <- function(u,phi,kappa) {
      u <- u+10e-16
      if(kappa==0.5) {
         out <- (u*(u-2*phi)*exp(-u/phi))/phi^4
      } else {
         bk <- besselK(u/phi,kappa)
         bk.p1 <- besselK(u/phi,kappa+1)
         bk.p2 <- besselK(u/phi,kappa+2)
         bk.m1 <- besselK(u/phi,kappa-1)
         bk.m2 <- besselK(u/phi,kappa-2)
         out <- (2^(-kappa-1)*phi^(-kappa-4)*u^kappa*(bk.p2*u^2+2*bk*u^2+
                bk.m2*u^2-4*kappa*bk.p1*phi*u-4*
                bk.p1*phi*u-4*kappa*bk.m1*phi*u-4*bk.m1*phi*u+
                4*kappa^2*bk*phi^2+4*kappa*bk*phi^2))/(gamma(kappa))
      }
      out
   }
    
   matern.grad.phi <- function(U,phi,kappa) {
    	n <- attr(U,"Size")
      grad.phi.mat <- matrix(NA,nrow=n,ncol=n)
      ind <- lower.tri(grad.phi.mat)             
      grad.phi <- der.phi(as.numeric(U),phi,kappa)
      grad.phi.mat[ind] <-  grad.phi
      grad.phi.mat <- t(grad.phi.mat)
      grad.phi.mat[ind] <-  grad.phi
      diag(grad.phi.mat) <- rep(der.phi(0,phi,kappa),n)
      grad.phi.mat
   }

   matern.hessian.phi <- function(U,phi,kappa) {
     	n <- attr(U,"Size")
	hess.phi.mat <- matrix(NA,nrow=n,ncol=n)
      ind <- lower.tri(hess.phi.mat)             
      hess.phi <- der2.phi(as.numeric(U),phi,kappa)
      hess.phi.mat[ind] <-  hess.phi
      hess.phi.mat <- t(hess.phi.mat)
      hess.phi.mat[ind] <-  hess.phi
      diag(hess.phi.mat) <- rep(der2.phi(0,phi,kappa),n)
      hess.phi.mat
   }
   
   if(length(ID.coords) > 0) {
      
      if(profile.llik) {
         profile.log.lik <- function(par) {
            phi <- par[1]
                  
            if(length(fixed.rel.nugget)==1) {
               nu2.1 <- fixed.rel.nugget
               nu2.2 <- par[2]
            } else {
               nu2.1 <- par[2]
               nu2.2 <- par[3]
            }
            R <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),
                 kappa=kappa)$varcov
            diag(R) <- diag(R)+nu2.1
            R.inv <- solve(R)
            Omega <- R.inv 
            diag(Omega) <- (diag(Omega)+n.coords/nu2.2)

            Omega.inv <- solve(Omega)
            M.beta <- DtD/nu2.2+
                -t(D.tilde)%*%Omega.inv%*%D.tilde/(nu2.2^2)
            v.beta <- Dty/nu2.2-t(D.tilde)%*%Omega.inv%*%y.tilde/(nu2.2^2)
            beta.hat <- as.numeric(solve(M.beta)%*%v.beta)
      
            M.det <- t(R*n.coords/nu2.2)
            diag(M.det) <- diag(M.det)+1
      
            mu.hat <- as.numeric(D%*%beta.hat)
            diff.y <- y-mu.hat 
            diff.y.tilde <- tapply(diff.y,ID.coords,sum)
            sigma2.hat <- as.numeric((sum(diff.y^2)/nu2.2-
                    t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde/(nu2.2^2))/
                    m)
     
            out <- -0.5*(m*log(sigma2.hat)+m*log(nu2.2)+
            as.numeric(determinant(M.det)$modulus))
            return(out)
         }
      
         compute.mle.std <- function(est.profile) {
            phi <- exp(est.profile$par[1])
            if(length(fixed.rel.nugget)==1) {
               nu2.1 <- fixed.rel.nugget
               nu2.2 <- exp(est.profile$par[2])
            } else {
               nu2.1 <- exp(est.profile$par[2])
               nu2.2 <- exp(est.profile$par[3])
            }
            R <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),
                 kappa=kappa)$varcov
            diag(R) <- diag(R)+nu2.1
            R.inv <- solve(R)
            Omega.star <- R.inv 
            diag(Omega.star) <- (diag(Omega.star)+n.coords/nu2.2)
            Omega.star.inv <- solve(Omega.star)
      
            M.beta <- DtD/nu2.2+
                   -t(D.tilde)%*%Omega.star.inv%*%D.tilde/(nu2.2^2)
            v.beta <- Dty/nu2.2-t(D.tilde)%*%Omega.star.inv%*%y.tilde/(nu2.2^2)
            beta <- as.numeric(solve(M.beta)%*%v.beta)
            
            mu <- as.numeric(D%*%beta)   
            diff.y <- y-mu  
            diff.y.tilde <- as.numeric(tapply(diff.y,ID.coords,sum))           
       
            v <- as.numeric(Omega.star.inv%*%diff.y.tilde)
            q.f1 <- sum(diff.y^2)
            q.f2 <- as.numeric(t(diff.y.tilde)%*%v)
            q.f <- q.f1/nu2.2-q.f2/(nu2.2^2)

            sigma2 <- as.numeric(q.f/m)

            grad.beta <- as.numeric(t(D)%*%(diff.y/nu2.2-
                           v[ID.coords]/(nu2.2^2)))/sigma2
            grad.sigma2 <- -0.5*(m-q.f/sigma2)

            M.det <- t(R*n.coords/nu2.2)
            diag(M.det) <- diag(M.det)+1
            M.det.inv <- solve(M.det)

            R1.phi <- matern.grad.phi(U,phi,kappa)*phi
            Omega.star.R.inv <- Omega.star.inv%*%R.inv
            der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
            der.M.phi <- t(R1.phi*n.coords/nu2.2)
            M1.trace.phi <- M.det.inv*t(der.M.phi)
            trace1.phi <- sum(M1.trace.phi)
            v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde)
            q.f.phi <-  as.numeric(t(diff.y.tilde)%*%v.beta.phi)
            grad.phi <- -0.5*(trace1.phi-
                  q.f.phi/((nu2.2^2)*sigma2))
   
            if(length(fixed.rel.nugget)==0) {
               der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
               trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords*nu2.1/nu2.2)
               v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                 diff.y.tilde)
               q.f.nu2.1 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.1)
               grad.nu2.1 <- -0.5*(trace1.nu2.1+
                  q.f.nu2.1/((nu2.2^2)*sigma2))
            }            
      
            der.Omega.star.inv.nu2.2 <- -t(Omega.star.inv*n.coords/nu2.2)%*%
                                       Omega.star.inv
            der.M.nu2.2 <- -t(R*n.coords/nu2.2)
            M1.trace.nu2.2 <- M.det.inv*t(der.M.nu2.2)
            trace1.nu2.2 <- sum(M1.trace.nu2.2)
            v.beta.nu2.2 <- as.numeric(der.Omega.star.inv.nu2.2%*%
                                 diff.y.tilde)
            q.f.nu2.2 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.2)
            grad.nu2.2 <- -0.5*(-q.f1/nu2.2+2*q.f2/(nu2.2^2))/sigma2+
                          -0.5*(m+trace1.nu2.2+
                                q.f.nu2.2/((nu2.2^2)*sigma2))

            ind.beta <- 1:p
            ind.sigma2 <- p+1
            ind.phi <- p+2
            if(length(fixed.rel.nugget)==0) {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.nu2.2)
               H <- matrix(NA,p+4,p+4)
               ind.nu2.1 <- p+3
               ind.nu2.2 <- p+4
            } else {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.2)
               H <- matrix(NA,p+3,p+3)
               ind.nu2.2 <- p+3
            }
            
      
            
            H[ind.beta,ind.beta] <- (-DtD/nu2.2+
                          t(D.tilde)%*%Omega.star.inv%*%D.tilde/(nu2.2^2))/sigma2

            H[ind.beta,ind.sigma2] <- H[ind.sigma2,ind.beta] <- -grad.beta

            H[ind.beta,ind.phi] <- H[ind.phi,ind.beta] <- -t(D)%*%v.beta.phi[ID.coords]/((nu2.2^2)*sigma2)      

            if(length(fixed.rel.nugget)==0) {

               H[ind.beta,ind.nu2.1] <- H[ind.nu2.1,ind.beta] <- t(D)%*%v.beta.nu2.1[ID.coords]/((nu2.2^2)*sigma2)      
 
            }
            
            H[ind.beta,ind.nu2.2] <- 
            H[ind.nu2.2,ind.beta] <- as.numeric(t(D)%*%(-diff.y/nu2.2+2*v[ID.coords]/(nu2.2^2))/sigma2)+
                                  t(D)%*%v.beta.nu2.2[ID.coords]/((nu2.2^2)*sigma2)

            H[ind.sigma2,ind.sigma2] <- -0.5*q.f/sigma2
            H[ind.sigma2,ind.phi] <- H[ind.phi,ind.sigma2] <- -(grad.phi+0.5*(trace1.phi))
            if(length(fixed.rel.nugget)==0) {
               H[ind.sigma2,ind.nu2.1] <- H[ind.nu2.1,ind.sigma2] <- -(grad.nu2.1+0.5*(trace1.nu2.1))
            }
            H[ind.sigma2,ind.nu2.2] <- H[ind.nu2.2,ind.sigma2] <- -(grad.nu2.2+0.5*(trace1.nu2.2+m))

            R2.phi <- matern.hessian.phi(U,phi,kappa)*(phi^2)+
                      R1.phi
            V1.phi <- -R.inv%*%R1.phi%*%R.inv
            V2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv
            der2.Omega.star.inv.phi <- Omega.star.inv%*%(2*V1.phi%*%
                                 Omega.star.inv%*%V1.phi-V2.phi)%*%
                                 Omega.star.inv
            M2.trace.phi <- M.det.inv*((R2.phi*n.coords/nu2.2))
            B2.phi <- M.det.inv%*%der.M.phi
      
            trace2.phi <- sum(M2.trace.phi)-
                    sum(B2.phi*t(B2.phi))
            H[ind.phi,ind.phi] <- -0.5*(trace2.phi-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi%*%diff.y.tilde)/((nu2.2^2)*sigma2))
   
            if(length(fixed.rel.nugget)==0) {

               V1.nu2.1 <- -R.inv%*%R.inv*nu2.1
               V2.phi.nu2.1 <- nu2.1*R.inv%*%(R1.phi%*%R.inv+R.inv%*%R1.phi)%*%R.inv
               der2.Omega.star.inv.phi.nu2.1 <- Omega.star.inv%*%(
                                       V1.phi%*%Omega.star.inv%*%V1.nu2.1+
                                       V1.nu2.1%*%Omega.star.inv%*%V1.phi-
                                       V2.phi.nu2.1)%*%
                                       Omega.star.inv
               B2.nu2.1 <- t(t(M.det.inv)*n.coords*nu2.1/nu2.2)
               trace2.phi.nu2.1 <- -sum(B2.phi*t(B2.nu2.1))

               H[ind.phi,ind.nu2.1] <- H[ind.nu2.1,ind.phi] <- -0.5*(trace2.phi.nu2.1-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))
   
            }
     
            der2.Omega.star.inv.phi.nu2.2 <- Omega.star.inv%*%(
                                       V1.phi%*%t(-Omega.star.inv*
                                       n.coords/nu2.2)+
                                       (-Omega.star.inv*
                                       n.coords/nu2.2)%*%V1.phi)%*%
                                       Omega.star.inv
            B2.nu2.2 <- M.det.inv%*%der.M.nu2.2
            trace2.phi.nu2.2 <- -trace1.phi-sum(B2.phi*t(B2.nu2.2))
            H[ind.phi,ind.nu2.2] <- H[ind.nu2.2,ind.phi] <- -0.5*(trace2.phi.nu2.2+
                    -as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                    +2*q.f.phi/((nu2.2^2)*sigma2))

            if(length(fixed.rel.nugget)==0) {
               aux.nu2.1 <- 2*R.inv*(nu2.1^2)
               diag(aux.nu2.1) <- diag(aux.nu2.1)-nu2.1
               V2.nu2.1 <- R.inv%*%aux.nu2.1%*%R.inv
               der2.Omega.star.inv.nu2.1 <- Omega.star.inv%*%(2*V1.nu2.1%*%
                                 Omega.star.inv%*%V1.nu2.1-V2.nu2.1)%*%
                                 Omega.star.inv
           
               trace2.nu2.1 <- trace1.nu2.1-
                    sum(B2.nu2.1*t(B2.nu2.1))

               H[ind.nu2.1,ind.nu2.1] <- -0.5*(trace2.nu2.1-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))
     
               der2.Omega.star.inv.nu2.1.nu2.2 <- Omega.star.inv%*%(
                                       V1.nu2.1%*%t(-Omega.star.inv*
                                       n.coords/nu2.2)+
                                       (-Omega.star.inv*
                                       n.coords/nu2.2)%*%V1.nu2.1)%*%
                                       Omega.star.inv
 
               trace2.nu2.1.nu2.2 <- -trace1.nu2.1-sum(B2.nu2.1*t(B2.nu2.2))
               H[ind.nu2.1,ind.nu2.2] <- H[ind.nu2.2,ind.nu2.1] <- -0.5*(trace2.nu2.1.nu2.2+
                    -as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.1.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                    -2*q.f.nu2.1/((nu2.2^2)*sigma2))

            }
         
            aux.nu2.2 <- 2*t(Omega.star.inv*n.coords)*n.coords/
                     (nu2.2^2)
            diag(aux.nu2.2) <- diag(aux.nu2.2)-n.coords/nu2.2
            der2.Omega.star.inv.nu2.2 <- -Omega.star.inv%*%(
                                   aux.nu2.2)%*%
                                   Omega.star.inv
      
            trace2.nu2.2 <- -trace1.nu2.2-
                    sum(B2.nu2.2*t(B2.nu2.2))
            H[ind.nu2.2,ind.nu2.2] <- -0.5*(q.f1/nu2.2-4*q.f2/(nu2.2^2)-4*q.f.nu2.2/(nu2.2^2))/sigma2+                 
                    -0.5*(trace2.nu2.2+as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2))
                    
      
            out <- list()
            if(length(fixed.rel.nugget)==0) {
               out$estimates <- c(beta,log(sigma2),log(phi),log(nu2.1),log(nu2.2))
            } else {
               out$estimates <- c(beta,log(sigma2),log(phi),log(nu2.2))
            }
            out$log.lik <- -est.profile$objective
            out$gradient <- g
            out$hessian <- H
            return(out)
         }

         estim <- list()           
         if(method=="nlminb") {
            est.profile <- nlminb(start.par,
               function(x) -profile.log.lik(exp(x)),
               control=list(trace=1*messages))
            est.summary <- compute.mle.std(est.profile)
            estim$estimate <- est.summary$estimates       
            if(messages) cat("\n","Gradient at MLE:",paste(round(est.summary$gradient,10)),"\n")
            estim$covariance <- solve(-est.summary$hessian)
            estim$log.lik <- est.summary$log.lik
         } else if(method=="BFGS") {
            est.profile <- maxBFGS(function(x) profile.log.lik(exp(x)),
            start=start.par,print.level=1*messages)
            est.profile$par <- est.profile$estimate
            est.profile$objective <- -est.profile$maximum
            est.summary <- compute.mle.std(est.profile)
            estim$estimate <- est.summary$estimates       
            if(messages) cat("Gradient at MLE: ",round(est.summary$gradient,10),"\n",sep=" ")
            estim$covariance <- solve(-est.summary$hessian)
            estim$log.lik <- est.summary$log.lik
         }
      } else {
         log.lik <- function(par) {
            beta <- par[1:p]
            sigma2 <- exp(par[p+1])
            phi <- exp(par[p+2])
            if(length(fixed.rel.nugget)==1) {
               nu2.1 <- fixed.rel.nugget
               nu2.2 <- exp(par[p+3])
            } else {
               nu2.1 <- exp(par[p+3])
               nu2.2 <- exp(par[p+4])
            }

            mu <- as.numeric(D%*%beta)   
            diff.y <- y-mu  
            diff.y.tilde <- as.numeric(tapply(diff.y,ID.coords,sum))

            R <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),
                 kappa=kappa)$varcov
            diag(R) <- diag(R)+nu2.1
            R.inv <- solve(R)
            Omega <- R.inv 
            diag(Omega) <- (diag(Omega)+n.coords/nu2.2)
            Omega <- Omega*(nu2.2^2)
            Omega.inv <- solve(Omega)
      
            q.f <- as.numeric(sum(diff.y^2)/nu2.2-
                   t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde)
            M.det <- t(R*n.coords/nu2.2)
            diag(M.det) <- diag(M.det)+1
            log.det.tot <- m*log(sigma2)+m*log(nu2.2)+
                     as.numeric(determinant(M.det)$modulus)
            out <- -0.5*(log.det.tot+q.f/sigma2)
            return(out)
         }
      
         grad.log.lik <- function(par) {
            beta <- par[1:p]
            sigma2 <- exp(par[p+1])
            phi <- exp(par[p+2])
            if(length(fixed.rel.nugget)==1) {
               nu2.1 <- fixed.rel.nugget
               nu2.2 <- exp(par[p+3])
            } else {
               nu2.1 <- exp(par[p+3])
               nu2.2 <- exp(par[p+4])
            }
            R <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),
                 kappa=kappa)$varcov
            diag(R) <- diag(R)+nu2.1
            R.inv <- solve(R)
            Omega.star <- R.inv 
            diag(Omega.star) <- (diag(Omega.star)+n.coords/nu2.2)
            Omega.star.inv <- solve(Omega.star)
      
            
            mu <- as.numeric(D%*%beta)   
            diff.y <- y-mu  
            diff.y.tilde <- as.numeric(tapply(diff.y,ID.coords,sum))           
       
            v <- as.numeric(Omega.star.inv%*%diff.y.tilde)
            q.f1 <- sum(diff.y^2)
            q.f2 <- as.numeric(t(diff.y.tilde)%*%v)
            q.f <- q.f1/nu2.2-q.f2/(nu2.2^2)
         
            grad.beta <- as.numeric(t(D)%*%(diff.y/nu2.2-
                           v[ID.coords]/(nu2.2^2)))/sigma2
            grad.sigma2 <- -0.5*(m-q.f/sigma2)

            M.det <- t(R*n.coords/nu2.2)
            diag(M.det) <- diag(M.det)+1
            M.det.inv <- solve(M.det)

            R1.phi <- matern.grad.phi(U,phi,kappa)*phi
            Omega.star.R.inv <- Omega.star.inv%*%R.inv
            der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
            der.M.phi <- t(R1.phi*n.coords/nu2.2)
            M1.trace.phi <- M.det.inv*t(der.M.phi)
            trace1.phi <- sum(M1.trace.phi)
            v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde)
            q.f.phi <-  as.numeric(t(diff.y.tilde)%*%v.beta.phi)
            grad.phi <- -0.5*(trace1.phi-
                  q.f.phi/((nu2.2^2)*sigma2))
   
            if(length(fixed.rel.nugget)==0) {
               der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
               trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords*nu2.1/nu2.2)
               v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                 diff.y.tilde)
               q.f.nu2.1 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.1)
               grad.nu2.1 <- -0.5*(trace1.nu2.1+
                  q.f.nu2.1/((nu2.2^2)*sigma2))
            }            
      
            der.Omega.star.inv.nu2.2 <- -t(Omega.star.inv*n.coords/nu2.2)%*%
                                       Omega.star.inv
            der.M.nu2.2 <- -t(R*n.coords/nu2.2)
            M1.trace.nu2.2 <- M.det.inv*t(der.M.nu2.2)
            trace1.nu2.2 <- sum(M1.trace.nu2.2)
            v.beta.nu2.2 <- as.numeric(der.Omega.star.inv.nu2.2%*%
                                 diff.y.tilde)
            q.f.nu2.2 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.2)
            grad.nu2.2 <- -0.5*(-q.f1/nu2.2+2*q.f2/(nu2.2^2))/sigma2+
                          -0.5*(m+trace1.nu2.2+
                                q.f.nu2.2/((nu2.2^2)*sigma2))

            ind.beta <- 1:p
            ind.sigma2 <- p+1
            ind.phi <- p+2
            if(length(fixed.rel.nugget)==0) {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.nu2.2)
            } else {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.2)
            }

            return(g)
         }
         
         hessian.log.lik <- function(par) {
            beta <- par[1:p]
            sigma2 <- exp(par[p+1])
            phi <- exp(par[p+2])
            if(length(fixed.rel.nugget)==1) {
               nu2.1 <- fixed.rel.nugget
               nu2.2 <- exp(par[p+3])
            } else {
               nu2.1 <- exp(par[p+3])
               nu2.2 <- exp(par[p+4])
            }
            R <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),
                 kappa=kappa)$varcov
            diag(R) <- diag(R)+nu2.1
            R.inv <- solve(R)
            Omega.star <- R.inv 
            diag(Omega.star) <- (diag(Omega.star)+n.coords/nu2.2)
            Omega.star.inv <- solve(Omega.star)
      
            mu <- as.numeric(D%*%beta)   
            diff.y <- y-mu  
            diff.y.tilde <- as.numeric(tapply(diff.y,ID.coords,sum))           
       
            v <- as.numeric(Omega.star.inv%*%diff.y.tilde)
            q.f1 <- sum(diff.y^2)
            q.f2 <- as.numeric(t(diff.y.tilde)%*%v)
            q.f <- q.f1/nu2.2-q.f2/(nu2.2^2)
         
            grad.beta <- as.numeric(t(D)%*%(diff.y/nu2.2-
                           v[ID.coords]/(nu2.2^2)))/sigma2
            grad.sigma2 <- -0.5*(m-q.f/sigma2)

            M.det <- t(R*n.coords/nu2.2)
            diag(M.det) <- diag(M.det)+1
            M.det.inv <- solve(M.det)

            R1.phi <- matern.grad.phi(U,phi,kappa)*phi
            Omega.star.R.inv <- Omega.star.inv%*%R.inv
            der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
            der.M.phi <- t(R1.phi*n.coords/nu2.2)
            M1.trace.phi <- M.det.inv*t(der.M.phi)
            trace1.phi <- sum(M1.trace.phi)
            v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde)
            q.f.phi <-  as.numeric(t(diff.y.tilde)%*%v.beta.phi)
            grad.phi <- -0.5*(trace1.phi-
                  q.f.phi/((nu2.2^2)*sigma2))
   
            if(length(fixed.rel.nugget)==0) {
               der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
               trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords*nu2.1/nu2.2)
               v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                 diff.y.tilde)
               q.f.nu2.1 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.1)
               grad.nu2.1 <- -0.5*(trace1.nu2.1+
                  q.f.nu2.1/((nu2.2^2)*sigma2))
            }            
      
            der.Omega.star.inv.nu2.2 <- -t(Omega.star.inv*n.coords/nu2.2)%*%
                                       Omega.star.inv
            der.M.nu2.2 <- -t(R*n.coords/nu2.2)
            M1.trace.nu2.2 <- M.det.inv*t(der.M.nu2.2)
            trace1.nu2.2 <- sum(M1.trace.nu2.2)
            v.beta.nu2.2 <- as.numeric(der.Omega.star.inv.nu2.2%*%
                                 diff.y.tilde)
            q.f.nu2.2 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.2)
            grad.nu2.2 <- -0.5*(-q.f1/nu2.2+2*q.f2/(nu2.2^2))/sigma2+
                          -0.5*(m+trace1.nu2.2+
                                q.f.nu2.2/((nu2.2^2)*sigma2))

            ind.beta <- 1:p
            ind.sigma2 <- p+1
            ind.phi <- p+2
            if(length(fixed.rel.nugget)==0) {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.nu2.2)
               H <- matrix(NA,p+4,p+4)
               ind.nu2.1 <- p+3
               ind.nu2.2 <- p+4
            } else {
               g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.2)
               H <- matrix(NA,p+3,p+3)
               ind.nu2.2 <- p+3
            }
            
      
            
            H[ind.beta,ind.beta] <- (-DtD/nu2.2+
                          t(D.tilde)%*%Omega.star.inv%*%D.tilde/(nu2.2^2))/sigma2

            H[ind.beta,ind.sigma2] <- H[ind.sigma2,ind.beta] <- -grad.beta

            H[ind.beta,ind.phi] <- H[ind.phi,ind.beta] <- -t(D)%*%v.beta.phi[ID.coords]/((nu2.2^2)*sigma2)      

            if(length(fixed.rel.nugget)==0) {

               H[ind.beta,ind.nu2.1] <- H[ind.nu2.1,ind.beta] <- t(D)%*%v.beta.nu2.1[ID.coords]/((nu2.2^2)*sigma2)      
 
            }
            
            H[ind.beta,ind.nu2.2] <- 
            H[ind.nu2.2,ind.beta] <- as.numeric(t(D)%*%(-diff.y/nu2.2+2*v[ID.coords]/(nu2.2^2))/sigma2)+
                                  t(D)%*%v.beta.nu2.2[ID.coords]/((nu2.2^2)*sigma2)

            H[ind.sigma2,ind.sigma2] <- -0.5*q.f/sigma2
            H[ind.sigma2,ind.phi] <- H[ind.phi,ind.sigma2] <- -(grad.phi+0.5*(trace1.phi))
            if(length(fixed.rel.nugget)==0) {
               H[ind.sigma2,ind.nu2.1] <- H[ind.nu2.1,ind.sigma2] <- -(grad.nu2.1+0.5*(trace1.nu2.1))
            }
            H[ind.sigma2,ind.nu2.2] <- H[ind.nu2.2,ind.sigma2] <- -(grad.nu2.2+0.5*(trace1.nu2.2+m))

            R2.phi <- matern.hessian.phi(U,phi,kappa)*(phi^2)+
                      R1.phi
            V1.phi <- -R.inv%*%R1.phi%*%R.inv
            V2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv
            der2.Omega.star.inv.phi <- Omega.star.inv%*%(2*V1.phi%*%
                                 Omega.star.inv%*%V1.phi-V2.phi)%*%
                                 Omega.star.inv
            M2.trace.phi <- M.det.inv*((R2.phi*n.coords/nu2.2))
            B2.phi <- M.det.inv%*%der.M.phi
      
            trace2.phi <- sum(M2.trace.phi)-
                    sum(B2.phi*t(B2.phi))
            H[ind.phi,ind.phi] <- -0.5*(trace2.phi-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi%*%diff.y.tilde)/((nu2.2^2)*sigma2))
   
            if(length(fixed.rel.nugget)==0) {

               V1.nu2.1 <- -R.inv%*%R.inv*nu2.1
               V2.phi.nu2.1 <- nu2.1*R.inv%*%(R1.phi%*%R.inv+R.inv%*%R1.phi)%*%R.inv
               der2.Omega.star.inv.phi.nu2.1 <- Omega.star.inv%*%(
                                       V1.phi%*%Omega.star.inv%*%V1.nu2.1+
                                       V1.nu2.1%*%Omega.star.inv%*%V1.phi-
                                       V2.phi.nu2.1)%*%
                                       Omega.star.inv
               B2.nu2.1 <- t(t(M.det.inv)*n.coords*nu2.1/nu2.2)
               trace2.phi.nu2.1 <- -sum(B2.phi*t(B2.nu2.1))

               H[ind.phi,ind.nu2.1] <- H[ind.nu2.1,ind.phi] <- -0.5*(trace2.phi.nu2.1-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))
   
            }
     
            der2.Omega.star.inv.phi.nu2.2 <- Omega.star.inv%*%(
                                       V1.phi%*%t(-Omega.star.inv*
                                       n.coords/nu2.2)+
                                       (-Omega.star.inv*
                                       n.coords/nu2.2)%*%V1.phi)%*%
                                       Omega.star.inv
            B2.nu2.2 <- M.det.inv%*%der.M.nu2.2
            trace2.phi.nu2.2 <- -trace1.phi-sum(B2.phi*t(B2.nu2.2))
            H[ind.phi,ind.nu2.2] <- H[ind.nu2.2,ind.phi] <- -0.5*(trace2.phi.nu2.2+
                    -as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.phi.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                    +2*q.f.phi/((nu2.2^2)*sigma2))

            if(length(fixed.rel.nugget)==0) {
               aux.nu2.1 <- 2*R.inv*(nu2.1^2)
               diag(aux.nu2.1) <- diag(aux.nu2.1)-nu2.1
               V2.nu2.1 <- R.inv%*%aux.nu2.1%*%R.inv
               der2.Omega.star.inv.nu2.1 <- Omega.star.inv%*%(2*V1.nu2.1%*%
                                 Omega.star.inv%*%V1.nu2.1-V2.nu2.1)%*%
                                 Omega.star.inv
           
               trace2.nu2.1 <- trace1.nu2.1-
                    sum(B2.nu2.1*t(B2.nu2.1))

               H[ind.nu2.1,ind.nu2.1] <- -0.5*(trace2.nu2.1-
                    as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))
     
               der2.Omega.star.inv.nu2.1.nu2.2 <- Omega.star.inv%*%(
                                       V1.nu2.1%*%t(-Omega.star.inv*
                                       n.coords/nu2.2)+
                                       (-Omega.star.inv*
                                       n.coords/nu2.2)%*%V1.nu2.1)%*%
                                       Omega.star.inv
 
               trace2.nu2.1.nu2.2 <- -trace1.nu2.1-sum(B2.nu2.1*t(B2.nu2.2))
               H[ind.nu2.1,ind.nu2.2] <- H[ind.nu2.2,ind.nu2.1] <- -0.5*(trace2.nu2.1.nu2.2+
                    -as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.1.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                    -2*q.f.nu2.1/((nu2.2^2)*sigma2))

            }
         
            aux.nu2.2 <- 2*t(Omega.star.inv*n.coords)*n.coords/
                     (nu2.2^2)
            diag(aux.nu2.2) <- diag(aux.nu2.2)-n.coords/nu2.2
            der2.Omega.star.inv.nu2.2 <- -Omega.star.inv%*%(
                                   aux.nu2.2)%*%
                                   Omega.star.inv
      
            trace2.nu2.2 <- -trace1.nu2.2-
                    sum(B2.nu2.2*t(B2.nu2.2))
            H[ind.nu2.2,ind.nu2.2] <- -0.5*(q.f1/nu2.2-4*q.f2/(nu2.2^2)-4*q.f.nu2.2/(nu2.2^2))/sigma2+                 
                    -0.5*(trace2.nu2.2+as.numeric(t(diff.y.tilde)%*%
                    der2.Omega.star.inv.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2))
                    
      
            return(H)
         }
         estim <- list()
         if(method=="nlminb") {
            estNLMINB <- nlminb(start.par,
                         function(x) -log.lik(x),
                         function(x) -grad.log.lik(x),
                         function(x) -hessian.log.lik(x),
                         control=list(trace=1*messages))

            estim$estimate <- estNLMINB$par
            hess.MLE <- hessian.log.lik(estim$estimate)
            if(messages) {
               grad.MLE <- grad.log.lik(estim$estimate)
               cat("\n","Gradient at MLE:",
                         paste(round(grad.MLE,10)),"\n")
            }
            estim$covariance <- solve(-hess.MLE)
            estim$log.lik <- -estNLMINB$objective

         } else if(method=="BFGS") {
            estimBFGS <- maxBFGS(log.lik,grad.log.lik,hessian.log.lik,
                         start.par,print.level=1*messages)   
            estim$estimate <- estimBFGS$estimate  
            if(messages) {
               cat("\n","Gradient at MLE:",
                         paste(round(estimBFGS$gradient,10)),"\n")
            }
            estim$covariance <- solve(-estimBFGS$hessian)
            estim$log.lik <- estimBFGS$maximum  
         }

      }         
   } else {
      log.lik <- function(par) {
	   beta <- par[1:p]
	   sigma2 <- exp(par[p+1])
	   phi <- exp(par[p+2])
	   if(length(fixed.rel.nugget)==1) {
            nu2 <- fixed.rel.nugget	
	   } else {
            nu2 <- exp(par[p+3])
         }
	   R <- varcov.spatial(dists.lowertri=U,kappa=kappa,
	                   cov.pars=c(1,phi),nugget=nu2)$varcov
         R.inv <- solve(R)
         ldet.R <- determinant(R)$modulus	                   
         mu <- as.numeric(D%*%beta)
         diff.y <- y-mu
	   out <- -0.5*(n*log(sigma2)+ldet.R+
	       t(diff.y)%*%R.inv%*%(diff.y)/sigma2)
	   as.numeric(out)
      }
	 
      grad.log.lik <- function(par) {
	   beta <- par[1:p]
         sigma2 <- exp(par[p+1])
         phi <- exp(par[p+2])
         if(length(fixed.rel.nugget)==1) {
            nu2 <- fixed.rel.nugget
	   } else {
            nu2 <- exp(par[p+3])
         }
		
	   R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                    cov.pars=c(1,phi),
	                    nugget=nu2,kappa=kappa)$varcov
	                                    
	   R.inv <- solve(R)
	   R1.phi <- matern.grad.phi(U,phi,kappa)                        
         m1.phi <- R.inv%*%R1.phi
         t1.phi <- -0.5*sum(diag(m1.phi))
         m2.phi <- m1.phi%*%R.inv; rm(m1.phi)           
         
         if(length(fixed.rel.nugget)==0) {
            t1.nu2 <- -0.5*sum(diag(R.inv))
            m2.nu2 <- R.inv%*%R.inv
         }               
         mu <- D%*%beta
         diff.y <- y-mu      
         q.f <- t(diff.y)%*%R.inv%*%diff.y
         grad.beta <-  t(D)%*%R.inv%*%(diff.y)/sigma2
         grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2 
         grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.y)%*%m2.phi%*%(diff.y))/sigma2)*phi
         if(length(fixed.rel.nugget)==0){
            grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/sigma2)*nu2
            out <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2)
         } else {
            out <- c(grad.beta,grad.log.sigma2,grad.log.phi)
         }        
         return(out)   
      }
	 
      hess.log.lik <- function(par) {
         beta <- par[1:p]
	   sigma2 <- exp(par[p+1])
	   phi <- exp(par[p+2])
	   mu <- D%*%beta
	   if(length(fixed.rel.nugget)==1) {
            nu2 <- fixed.rel.nugget	
	   } else {
	      nu2 <- exp(par[p+3])
	   }
		
   	   R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                    cov.pars=c(1,phi),
	                    nugget=nu2,kappa=kappa)$varcov
	                                    
	   R.inv <- solve(R)
	   R1.phi <- matern.grad.phi(U,phi,kappa)                        
         m1.phi <- R.inv%*%R1.phi
         t1.phi <- -0.5*sum(diag(m1.phi))
         m2.phi <- m1.phi%*%R.inv; rm(m1.phi)           
         
         if(length(fixed.rel.nugget)==0) {    
            t1.nu2 <- -0.5*sum(diag(R.inv))
            m2.nu2 <- R.inv%*%R.inv   
            t2.nu2 <- 0.5*sum(diag(m2.nu2))
            n2.nu2 <- 2*R.inv%*%m2.nu2
            t2.nu2.phi <- 0.5*sum(diag(R.inv%*%R1.phi%*%R.inv))
            n2.nu2.phi <- R.inv%*%(R.inv%*%R1.phi+
                        R1.phi%*%R.inv)%*%R.inv  
         }
         
         R2.phi <- matern.hessian.phi(U,phi,kappa)
         t2.phi <- -0.5*sum(diag(R.inv%*%R2.phi-R.inv%*%R1.phi%*%R.inv%*%R1.phi))
         n2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv          
                     
         diff.y <- y-mu      
         q.f <- t(diff.y)%*%R.inv%*%diff.y
         grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2 
         grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.y)%*%m2.phi%*%(diff.y))/sigma2)*phi
         if(length(fixed.rel.nugget)==0){
            grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/sigma2)*nu2
         }
         
         H <- matrix(0,nrow=length(par),ncol=length(par)) 
         H[1:p,1:p] <- -t(D)%*%R.inv%*%D/sigma2    
         H[1:p,p+1] <- H[p+1,1:p] <- -t(D)%*%R.inv%*%(diff.y)/sigma2
         H[1:p,p+2] <- H[p+2,1:p] <- -phi*as.numeric(t(D)%*%m2.phi%*%(diff.y))/sigma2
         H[p+1,p+1] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
                                grad.log.sigma2
                                
         H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
                                         
         H[p+2,p+2] <- (t2.phi-0.5*t(diff.y)%*%n2.phi%*%(diff.y)/sigma2)*phi^2+
                                            grad.log.phi          
                                            
         if(length(fixed.rel.nugget)==0) {
            H[1:p,p+3] <- H[p+3,1:p] <- -nu2*as.numeric(t(D)%*%m2.nu2%*%(diff.y))/sigma2 
            H[p+2,p+3] <- H[p+3,p+2] <- (t2.nu2.phi-0.5*t(diff.y)%*%n2.nu2.phi%*%(diff.y)/sigma2)*phi*nu2
            H[p+1,p+3] <- H[p+3,p+1] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
            H[p+3,p+3] <- (t2.nu2-0.5*t(diff.y)%*%n2.nu2%*%(diff.y)/sigma2)*nu2^2+
                                         grad.log.nu2
         }	                   
         return(H)               
      }
      estim <- list()     
      if(method=="BFGS") {
         estimBFGS <- maxBFGS(log.lik,grad.log.lik,hess.log.lik,
                           start.par,print.level=1*messages)   
         estim$estimate <- estimBFGS$estimate  
         estim$covariance <- solve(-estimBFGS$hessian)
         estim$log.lik <- estimBFGS$maximum                                                 
      } else if(method=="nlminb") {
         estimNLMINB <- nlminb(start.par,function(x) -log.lik(x),
         function(x) -grad.log.lik(x),
         function(x) -hess.log.lik(x),control=list(trace=1*messages)) 
         estim$estimate <- estimNLMINB$par            
         estim$covariance <- solve(-hess.log.lik(estimNLMINB$par))
         estim$log.lik <- -estimNLMINB$objective              
      }
   }
   names(estim$estimate)[1:p] <- colnames(D)
   names(estim$estimate)[(p+1):(p+2)] <- c("log(sigma^2)","log(phi)")
   if(length(fixed.rel.nugget)==0) names(estim$estimate)[p+3] <- "log(nu^2)"
   if(length(ID.coords)>0) {
      if(length(fixed.rel.nugget)==0) {
         names(estim$estimate)[p+4] <- "log(nu.star^2)"
      } else {
         names(estim$estimate)[p+3] <- "log(nu.star^2)"
      }
   }
   rownames(estim$covariance) <- colnames(estim$covariance) <- names(estim$estimate)    
   estim$y <- y
   estim$D <- D
   estim$coords <- coords
   if(length(ID.coords)>0) estim$ID.coords <- ID.coords
   estim$method <- method 
   estim$kappa <- kappa       
   if(length(fixed.rel.nugget)==1) {
      estim$fixed.rel.nugget <- fixed.rel.nugget
   } else {
      estim$fixed.rel.nugget <- NULL
   }	
   class(estim) <- "PrevMap"
   return(estim)
}

##' @title Spatial predictions for the geostatistical Linear Gaussian model using plug-in of ML estimates
##' @description This function performs spatial prediction, fixing the model parameters at the maximum likelihood estimates of a linear geostatistical model.
##' @param object an object of class "PrevMap" obtained as result of a call to \code{\link{linear.model.MLE}}.
##' @param grid.pred a matrix of prediction locations.
##' @param predictors a data frame of the values of the explanatory variables at each of the locations in \code{grid.pred}; each column correspond to a variable and each row to a location. \bold{Warning:} the names of the columns in the data frame must match those in the data used to fit the model. Default is \code{predictors=NULL} for models with only an intercept.
##' @param type a character indicating the type of spatial predictions: \code{type="marginal"} for marginal predictions or \code{type="joint"} for joint predictions. Default is \code{type="marginal"}. In the case of a low-rank approximation only marginal predictions are available.
##' @param scale.predictions a character vector of maximum length 3, indicating the required scale on which spatial prediction is carried out: "logit", "prevalence" and "odds". Default is \code{scale.predictions=c("logit","prevalence","odds")}.
##' @param quantiles a vector of quantiles used to summarise the spatial predictions.
##' @param n.sim.prev number of simulation for predictions of prevalence. Default is \code{n.sim.prev=1000}.
##' @param standard.errors logical; if \code{standard.errors=TRUE}, then standard errors for each \code{scale.predictions} are returned. Default is \code{standard.errors=FALSE}.
##' @param thresholds a vector of exceedance thresholds; default is \code{thresholds=NULL}.
##' @param scale.thresholds a character value indicating the scale on which exceedance thresholds are provided; \code{"logit"}, \code{"prevalence"} or \code{"odds"}. Default is \code{scale.thresholds=NULL}.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param include.nugget logical; if \code{include.nugget=TRUE} then the nugget effect is included in the predictions. This option is available only for fitted linear models with locations having multiple observations. Default is \code{include.nugget=FALSE}.
##' @return A "pred.PrevMap" object list with the following components: \code{logit}; \code{prevalence}; \code{odds}; \code{exceedance.prob}, corresponding to a matrix of the exceedance probabilities where each column corresponds to a specified value in \code{thresholds}; \code{grid.pred} prediction locations; \code{samples}, corresponding to the predictive samples of the linear predictor (only if \code{any(scale.predictions=="prevalence")}). 
##' Each of the three components \code{logit}, \code{prevalence} and  \code{odds} is also a list with the following components:
##' @return \code{predictions}: a vector of the predictive mean for the associated quantity (logit, odds or prevalence).
##' @return \code{standard.errors}: a vector of prediction standard errors (if \code{standard.errors=TRUE}).
##' @return \code{quantiles}: a matrix of quantiles of the resulting predictions with each column corresponding to a quantile specified through the argument \code{quantiles}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}  
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom geoR matern 
##' @export 
spatial.pred.linear.MLE <- function(object,grid.pred,predictors=NULL,
                                    type="marginal",
                                    scale.predictions=c("logit","prevalence","odds"),
                                    quantiles=c(0.025,0.975),n.sim.prev=1000,
                                    standard.errors=FALSE,
                                    thresholds=NULL,scale.thresholds=NULL,
                                    messages=TRUE,include.nugget=FALSE) {
   linear.ID.coords <- length(object$ID.coords) > 0 & 
                       substr(object$call[1],1,6)=="linear"
   if(linear.ID.coords & messages & !include.nugget) cat("NOTE: the nugget effect IS NOT included in the predictions. \n")
   if(linear.ID.coords & messages & include.nugget) cat("NOTE: the nugget effect IS included in the predictions. \n")
   if(include.nugget & !linear.ID.coords) warning("the argument 'include.nugget' is ignored; the option of
   including the nugget effect in the target predictions is 
   available only for models having locations with multiple observations.")   
   if(nrow(grid.pred) < 2) stop("prediction locations must be at least two.")
   if(length(predictors)>0 && class(predictors)!="data.frame") stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.")
   if(length(predictors)>0 && any(is.na(predictors))) stop("missing values found in 'predictors'.")
   p <- ncol(object$D)
   kappa <- object$kappa	
   n.pred <- nrow(grid.pred)
   coords <- object$coords
   out <- list()
   
   for(i in 1:length(scale.predictions)) {
      if(any(c("logit","prevalence","odds")==scale.predictions[i])==
		          FALSE) stop("invalid scale.predictions")
   }
   if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions must be marginal or joint.")

	
   if(length(thresholds)>0) {
      if(any(c("logit","prevalence","odds")==scale.thresholds)==FALSE) {
	   stop("scale thresholds must be logit, prevalence or odds scale.")
	}
   }
   if(length(thresholds)==0 & length(scale.thresholds)>0 |
      length(thresholds)>0 & length(scale.thresholds)==0) stop("to estimate exceedance probabilities both thresholds and scale.thresholds.")

   if(p==1) {
	predictors <- matrix(1,nrow=n.pred)	
   } else {
	if(length(dim(predictors))==0) stop("covariates at prediction locations should be provided.")	
	predictors <- as.matrix(model.matrix(delete.response(terms(formula(object$call))),data=predictors))
	if(nrow(predictors)!=nrow(grid.pred)) stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")
	if(ncol(predictors)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
   }
        
   beta <- object$estimate[1:p]
   mu <- as.numeric(object$D%*%beta)
   ck <- length(dim(object$knots)) > 0
   if(ck & type=="joint") {
      warning("only marginal predictions are available for the low-rank approximation.")
    	type<-"marginal"
   }
    
   if(ck) {	
      sigma2 <- exp(object$estimate[p+1])/object$const.sigma2
	rho <- 2*sqrt(object$kappa)*exp(object$estimate[p+2])
	if(length(object$fixed.rel.nugget)==0){		
	   tau2 <- sigma2*exp(object$estimate[p+3]) 
	} else {
	   tau2 <- object$fixed.rel.nugget*sigma2
	}
    
      inv.vcov <- function(nu2,A,K) {
	   mat <- (nu2^(-1))*A
	   diag(mat) <- diag(mat)+1
	   mat <- solve(mat)
	   mat <- -(nu2^(-1))*K%*%mat%*%t(K)
         diag(mat) <- diag(mat)+1
         mat*(nu2^(-1))
      }

      mu.pred <- as.numeric(predictors%*%beta)  
      U.k <- as.matrix(pdist(coords,object$knots))
	U.k.pred <- as.matrix(pdist(grid.pred,object$knots))
	K <- matern.kernel(U.k,rho,kappa)
	K.pred <- matern.kernel(U.k.pred,rho,kappa)
      
      C <- K.pred%*%t(K)*sigma2
      Sigma.inv <- inv.vcov(tau2/sigma2,t(K)%*%K,K)/sigma2
      A <- C%*%Sigma.inv
      
      mu.cond <- as.vector(mu.pred+C%*%Sigma.inv%*%(object$y-mu))
      sd.cond <- sqrt(sigma2*apply(K.pred,1,function(x) sum(x^2))-
                 apply(A*C,1,sum))
   } else {	
	sigma2 <- exp(object$estimate[p+1])
	phi <- exp(object$estimate[p+2])
	if(length(object$fixed.rel.nugget)==0){		
	   tau2 <- sigma2*exp(object$estimate[p+3]) 
         if(linear.ID.coords) omega2 <- sigma2*exp(object$estimate[p+4]) 
	} else {
	   tau2 <- object$fixed.rel.nugget*sigma2
         if(linear.ID.coords) omega2 <- sigma2*exp(object$estimate[p+3])
	}
      mu.pred <- as.numeric(predictors%*%beta)  

	U <- dist(coords)
	U.pred.coords <- as.matrix(pdist(grid.pred,coords))	
	Sigma <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov 
	Sigma.inv <- solve(Sigma)   	      
      C <- sigma2*matern(U.pred.coords,phi,kappa)	                
      if(linear.ID.coords) {
         nu2.2 <- omega2/sigma2
         n.coords <- as.numeric(tapply(object$ID.coords,object$ID.coords,length))
         diff.y <- object$y-mu
         diff.y.tilde <- tapply(diff.y, object$ID.coords,sum)
         Omega.star <- Sigma.inv*sigma2
         diag(Omega.star) <- diag(Omega.star)+n.coords/nu2.2
         Omega.star.inv <- solve(Omega.star)
         M.sd <- Omega.star.inv%*%(t(C)*n.coords)
         mu.cond <- mu.pred+C%*%diff.y.tilde/(sigma2*nu2.2)+
                    -t(t(C)*n.coords)%*%Omega.star.inv%*%diff.y.tilde/
                   (sigma2*nu2.2^2)

         if(type=="marginal") {
            if(include.nugget) {
               sd.cond <- sqrt(sigma2+tau2-
                       apply(
                       (t(t(C)*n.coords)*C)/(sigma2*nu2.2)+
                       -t(t(C)*n.coords)*t(M.sd)/(sigma2*nu2.2^2),1,sum))

            } else {
               sd.cond <- sqrt(sigma2-
                          apply(
                          (t(t(C)*n.coords)*C)/(sigma2*nu2.2)+
                          -t(t(C)*n.coords)*t(M.sd)/(sigma2*nu2.2^2),1,sum))
            } 
         }
      } else {
         A <- C%*%Sigma.inv
		       	             	    		 	   	
         mu.cond <- as.numeric(mu.pred+A%*%(object$y-mu))    
	   if(type=="marginal") sd.cond <- sqrt(sigma2-apply(A*C,1,sum))
      }
   } 

   if(type=="joint" & any(scale.predictions=="prevalence")) {
      if(messages) cat("Type of prevalence predictions: joint (this step might be demanding) \n")
      if(linear.ID.coords) {
         if(include.nugget) {
            Sigma.pred <-  varcov.spatial(coords=grid.pred,
                           cov.model="matern",nugget=tau2,
	                     cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
         } else {
            Sigma.pred <-  varcov.spatial(coords=grid.pred,
                           cov.model="matern",nugget=0,
	                     cov.pars=c(sigma2,phi),kappa=kappa)$varcov 

         }
         Sigma.cond <- Sigma.pred+
                       -t(t(C)*n.coords)%*%t(C)/(sigma2*nu2.2)+
                       t(t(C)*n.coords)%*%M.sd/(sigma2*nu2.2^2)
         sd.cond <- sqrt(diag(Sigma.cond))
      } else {
         A <- C%*%Sigma.inv
    	   Sigma.pred <-  varcov.spatial(coords=grid.pred,cov.model="matern",
	                cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
         Sigma.cond <- Sigma.pred - A%*%t(C)
	   sd.cond <- sqrt(diag(Sigma.cond))	   
      }
   }  
    
   if(any(scale.predictions=="prevalence")) {
      if(type=="marginal") {
         if(messages) cat("Type of prevalence predictions: marginal \n")
         eta.sim <- sapply(1:n.sim.prev, function(i) rnorm(n.pred,mu.cond,sd.cond))
      } else if(type=="joint") {
         Sigma.cond.sroot <- t(chol(Sigma.cond))
         eta.sim <- sapply(1:n.sim.prev, function(i) mu.cond+Sigma.cond.sroot%*%rnorm(n.pred))                                                                   
      }      	
   }
   
   if(any(scale.predictions=="logit")) {
      if(messages) cat("Spatial predictions: logit \n")    	   
    	out$logit$predictions <-  mu.cond
      if(standard.errors) {
         out$logit$standard.errors <- sd.cond
      }
      if(length(quantiles) > 0) {
         out$logit$quantiles <- sapply(quantiles,function(q) qnorm(q,mean=mu.cond,sd=sd.cond))
      }
      
      if(length(thresholds) > 0 && scale.thresholds=="logit") {
         out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
         out$exceedance.prob <- sapply(thresholds, function(x)
                                                           1-pnorm(x,mean=mu.cond,sd=sd.cond))	      
         colnames(out$exceedance.prob) <- paste(thresholds,sep="")                                                                                                                
      }		   
   }	
              
    
   if(any(scale.predictions=="prevalence")) {
      if(messages) cat("Spatial predictions: prevalence \n") 
      prev.sim <- 1/(1+exp(-eta.sim))  	   
      out$prevalence$predictions <- apply(prev.sim,1,mean)
      if(standard.errors) {
         out$prevalence$standard.errors <- apply(prev.sim,1,sd)
      }
       
      if(length(quantiles) > 0) {
         out$prevalence$quantiles <- t(apply(prev.sim,1,
                                     function(r) quantile(r,quantiles)))         
      }
       	
      if(length(thresholds) > 0 && scale.thresholds=="prevalence") {       	
         out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
         out$exceedance.prob <- sapply(log(thresholds/(1-thresholds)), function(x)
                                                           1-pnorm(x,mean=mu.cond,sd=sd.cond))	 
         colnames(out$exceedance.prob) <- paste(thresholds,sep="")                                                           
      }	   
   }
   
   if(any(scale.predictions=="odds")) {
      if(messages) cat("Spatial predictions: odds \n")           	   
      out$odds$predictions <- exp(mu.cond+0.5*(sd.cond^2))
      if(standard.errors) {
         out$odds$standard.errors <- sqrt(exp(2*mu.cond+sd.cond^2)*(exp(sd.cond^2)-1))
      }
       
      if(length(quantiles) > 0) {
         out$odds$quantiles <- sapply(quantiles,function(q) qlnorm(q,meanlog=mu.cond,sdlog=sd.cond))         
      }
       	
      if(length(thresholds) > 0 && scale.thresholds=="odds") {
         out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))       	  
         out$exceedance.prob <- sapply(log(thresholds), function(x)
                                                           1-pnorm(x,mean=mu.cond,sd=sd.cond))	 
         colnames(out$exceedance.prob) <- paste(thresholds,sep="")                                                           
      }	   
   }  
   out$grid.pred <- grid.pred
   if(any(scale.predictions=="prevalence")) {
      out$samples <- eta.sim
   }
   class(out) <- "pred.PrevMap"
   out
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
##' @importFrom pdist pdist
geo.linear.MLE.lr <- function(formula,coords,knots,data,kappa,start.cov.pars,
                                           method="BFGS",messages=TRUE) {          
    knots <- as.matrix(knots)
    start.cov.pars <- as.numeric(start.cov.pars)
    kappa <- as.numeric(kappa)
    	coords <- as.matrix(model.frame(coords,data))
	if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")    
	if(any(method==c("BFGS","nlminb"))==FALSE) stop("method must be either BFGS or nlminb.")                                       	                                     	
	mf <- model.frame(formula,data=data)
	y <- as.numeric(model.response(mf))
	n <- length(y)
	D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
	p <- ncol(D)

	if(length(start.cov.pars)!=2) stop("wrong length of start.cov.pars")		
		    
    U.k <- as.matrix(pdist(coords,knots))
    
    log.det.vcov <- function(nu2,A) {
	   mat <- A/nu2
	   diag(mat) <- diag(mat)+1
	   determinant(mat)$modulus+(n*log(nu2))
    }

    inv.vcov <- function(nu2,A,K) {
	   mat <- (nu2^(-1))*A
	   diag(mat) <- diag(mat)+1
	   mat <- solve(mat)
	   mat <- -(nu2^(-1))*K%*%mat%*%t(K)
       diag(mat) <- diag(mat)+1
       mat*(nu2^(-1))
    }
     
    if(any(start.cov.pars<0)) stop("start.cov.pars must be positive.")
    start.cov.pars[1] <- 2*sqrt(kappa)*start.cov.pars[1]
    K.start <- matern.kernel(U.k,start.cov.pars[1],kappa)
    A.start <- t(K.start)%*%K.start
    R.inv <- inv.vcov(start.cov.pars[2],A.start,K.start)
    beta.start <- as.numeric(solve(t(D)%*%R.inv%*%D)%*%t(D)%*%R.inv%*%y)
    diff.b <- y-D%*%beta.start
    sigma2.start <- as.numeric(t(diff.b)%*%R.inv%*%diff.b/length(y))
    start.par <- c(beta.start,sigma2.start,start.cov.pars) 
    
    start.par[-(1:p)] <- log(start.par[-(1:p)])
    der.rho <- function(u,rho,kappa) {
        u <- u+10e-16
        if(kappa==2) {
           out <- ((2^(7/2)*u-4*rho)*exp(-(2^(3/2)*u)/rho))/(sqrt(pi)*rho^3)
        } else {
           out <- (kappa^(kappa/4+1/4)*sqrt(gamma(kappa+1))*rho^(-kappa/2-5/2)*u^(kappa/2-1/2)*
                     (2*sqrt(kappa)*besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*
                     u+2*besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*sqrt(kappa)*
                     u-besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*rho-
                     besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*rho))/
                     (sqrt(gamma(kappa))*gamma((kappa+1)/2)*sqrt(pi))
        }
        out
     }

    der2.rho <- function(u,rho,kappa) {
        u <- u+10e-16
        if(kappa==2) {
          	out <- (8*(4*u^2-2^(5/2)*rho*u+rho^2)*exp(-(2^(3/2)*u)/rho))/(sqrt(pi)*rho^5)
        } else {
            out <- (kappa^(kappa/4+1/4)*sqrt(gamma(kappa+1))*rho^(-kappa/2-9/2)*u^(kappa/2-1/2)*
                      (4*kappa*besselK((2*sqrt(kappa)*u)/rho,(kappa+3)/2)*u^2+
                      8*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*u^2+
                      4*besselK((2*sqrt(kappa)*u)/rho,(kappa-5)/2)*kappa*u^2-
                      4*kappa^(3/2)*besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*rho*u-12*sqrt(kappa)*
                      besselK((2*sqrt(kappa)*u)/rho,(kappa+1)/2)*rho*u-4*
                      besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*kappa^(3/2)*rho*u-
                      12*besselK((2*sqrt(kappa)*u)/rho,(kappa-3)/2)*sqrt(kappa)*rho*u+
                      besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa^2*rho^2+
                      4*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*kappa*rho^2+
                      3*besselK((2*sqrt(kappa)*u)/rho,(kappa-1)/2)*rho^2))/
                     (2*sqrt(gamma(kappa))*gamma((kappa+1)/2)*sqrt(pi))
        }
        out
    }
         
	log.lik <- function(par) {
		beta <- par[1:p]
		sigma2 <- exp(par[p+1])
		rho <- exp(par[p+2])
		nu2 <- exp(par[p+3])
		
	    K <- matern.kernel(U.k,rho,kappa)
	    A <- t(K)%*%K
        R.inv <- inv.vcov(nu2,A,K)
        ldet.R <- log.det.vcov(nu2,A)	                   
        mu <- as.numeric(D%*%beta)
        diff.y <- y-mu
	    out <- -0.5*(n*log(sigma2)+ldet.R+
	               t(diff.y)%*%R.inv%*%(diff.y)/sigma2)
	    as.numeric(out)
	 }
	 
	 grad.log.lik <- function(par) {
	 	beta <- par[1:p]
		sigma2 <- exp(par[p+1])
		rho <- exp(par[p+2])
		nu2 <- 	exp(par[p+3])
		
		K <- matern.kernel(U.k,rho,kappa)
	    A <- t(K)%*%K                                
	    R.inv <- inv.vcov(nu2,A,K)
	    K.d <- der.rho(U.k,rho,kappa)
	    R1.rho <- K.d%*%t(K)+K%*%t(K.d)                        
        m1.rho <- R.inv%*%R1.rho
        t1.rho <- -0.5*sum(diag(m1.rho))
        m2.rho <- m1.rho%*%R.inv; rm(m1.rho)           
         

        t1.nu2 <- -0.5*sum(diag(R.inv))
        m2.nu2 <- R.inv%*%R.inv
        
        mu <- D%*%beta
        diff.y <- y-mu      
        q.f <- t(diff.y)%*%R.inv%*%diff.y
        grad.beta <-  t(D)%*%R.inv%*%(diff.y)/sigma2
        grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2 
        grad.log.rho <- (t1.rho+0.5*as.numeric(t(diff.y)%*%m2.rho%*%(diff.y))/sigma2)*rho
        grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/
                                    sigma2)*nu2
        out <- c(grad.beta,grad.log.sigma2,grad.log.rho,grad.log.nu2)    
        return(out)   
	 }
			
	 hess.log.lik <- function(par) {
	 	beta <- par[1:p]
		sigma2 <- exp(par[p+1])
		rho <- exp(par[p+2])
		nu2 <- 	exp(par[p+3])

		
		K <- matern.kernel(U.k,rho,kappa)
	    A <- t(K)%*%K                                
	    R.inv <- inv.vcov(nu2,A,K)
	    K.d <- der.rho(U.k,rho,kappa)
	    R1.rho <- K.d%*%t(K)+K%*%t(K.d)                        
        m1.rho <- R.inv%*%R1.rho
        t1.rho <- -0.5*sum(diag(m1.rho))
        m2.rho <- m1.rho%*%R.inv; rm(m1.rho)           
        
        mu <- D%*%beta
        diff.y <- y-mu    
         
        t1.nu2 <- -0.5*sum(diag(R.inv))
        m2.nu2 <- R.inv%*%R.inv
   
        t1.nu2 <- -0.5*sum(diag(R.inv))
        m2.nu2 <- R.inv%*%R.inv   
        t2.nu2 <- 0.5*sum(diag(m2.nu2))
        n2.nu2 <- 2*R.inv%*%m2.nu2
        t2.nu2.rho <- 0.5*sum(diag(R.inv%*%R1.rho%*%R.inv))
        n2.nu2.rho <- R.inv%*%(R.inv%*%R1.rho+
                               R1.rho%*%R.inv)%*%R.inv  
         
         K.d2 <- der2.rho(U.k,rho,kappa)
         R2.rho <- K.d2%*%t(K)+K%*%t(K.d2)+2*K.d%*%t(K.d)
         t2.rho <- -0.5*sum(diag(R.inv%*%R2.rho-R.inv%*%R1.rho%*%R.inv%*%R1.rho))
         n2.rho <- R.inv%*%(2*R1.rho%*%R.inv%*%R1.rho-R2.rho)%*%R.inv          
                       
         q.f <- t(diff.y)%*%R.inv%*%diff.y
         grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2 
         grad.log.rho <- (t1.rho+0.5*as.numeric(t(diff.y)%*%m2.rho%*%(diff.y))/
                                   sigma2)*rho        
          grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/
                                      sigma2)*nu2
         
         H <- matrix(0,nrow=length(par),ncol=length(par)) 
         H[1:p,1:p] <- -t(D)%*%R.inv%*%D/sigma2    
         H[1:p,p+1] <- H[p+1,1:p] <- -t(D)%*%R.inv%*%(diff.y)/sigma2
         H[1:p,p+2] <- H[p+2,1:p] <- -rho*as.numeric(t(D)%*%m2.rho%*%(diff.y))/sigma2
         H[p+1,p+1] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
                                grad.log.sigma2
                                
         H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.rho/rho-t1.rho)*(-rho)
                                         
         H[p+2,p+2] <- (t2.rho-0.5*t(diff.y)%*%n2.rho%*%(diff.y)/sigma2)*rho^2+
                                grad.log.rho          
                                            

         	H[1:p,p+3] <- H[p+3,1:p] <- -nu2*as.numeric(t(D)%*%m2.nu2%*%(diff.y))/
         	                                         sigma2 
         	H[p+2,p+3] <- H[p+3,p+2] <- (t2.nu2.rho-0.5*t(diff.y)%*%n2.nu2.rho%*%(diff.y)/
         	                                               sigma2)*rho*nu2
         	H[p+1,p+3] <- H[p+3,p+1] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
         	H[p+3,p+3] <- (t2.nu2-0.5*t(diff.y)%*%n2.nu2%*%(diff.y)/sigma2)*nu2^2+
                                         grad.log.nu2
                   
         return(H)               
	 }
	 	
     estim <- list()     
     if(method=="BFGS") {
        estimBFGS <- maxBFGS(log.lik,grad.log.lik,hess.log.lik,
                           start.par,print.level=1*messages)   
        estim$estimate <- estimBFGS$estimate  
        hessian <- estimBFGS$hessian
        estim$log.lik <- estimBFGS$maximum                                                 
     }
      
      if(method=="nlminb") {
      	estimNLMINB <- nlminb(start.par,function(x) -log.lik(x),
      	            function(x) -grad.log.lik(x),
      	            function(x) -hess.log.lik(x),control=list(trace=1*messages)) 
      	estim$estimate <- estimNLMINB$par            
      	hessian <- hess.log.lik(estimNLMINB$par)
      	estim$log.lik <- -estimNLMINB$objective              
      }                  
      K.hat <- matern.kernel(U.k,exp(estim$estimate[p+2]),kappa)
      const.sigma2 <- mean(apply(K.hat,1,function(r) sqrt(sum(r^2))))
      const.sigma2.grad <- mean(apply(U.k,1,function(r) {
      	                          rho <- exp(estim$estimate[p+2])
		                          k <- matern.kernel(r,rho,kappa)
		                          k.grad <- der.rho(r,rho,kappa)
		                          den <- sqrt(sum(k^2))
		                          num <- sum(k*k.grad)
	                              num/den
	                              }))
      estim$estimate[p+1] <- estim$estimate[p+1]+log(const.sigma2)
      estim$estimate[p+2] <- estim$estimate[p+2]-log(2*sqrt(kappa))      
      J <- diag(1,p+3)
      J[p+1,p+2] <- exp(estim$estimate[p+2])*const.sigma2.grad/const.sigma2
      hessian <- J%*%hessian%*%t(J) 
      names(estim$estimate)[1:p] <- colnames(D)
      names(estim$estimate)[(p+1):(p+3)] <- c("log(sigma^2)","log(phi)","log(nu^2)")
      estim$covariance <- solve(-J%*%hessian%*%t(J))
      colnames(estim$covariance) <- rownames(estim$covariance) <- 
      names(estim$estimate)     
      estim$y <- y
      estim$D <- D
      estim$coords <- coords
      estim$method <- method       
      estim$kappa <- kappa       
      estim$knots <- knots
      estim$fixed.rel.nugget <- NULL
      estim$const.sigma2 <- const.sigma2
      class(estim) <- "PrevMap"
      return(estim)
}


##' @title Priors specification 
##' @description This function is used to define priors for the model parameters of a Bayesian geostatistical model.
##' @param beta.mean mean vector of the Gaussian prior for the regression coefficients.
##' @param beta.covar covariance matrix of the Gaussian prior for the regression coefficients.
##' @param log.prior.sigma2 a function corresponding to the log-density of the prior distribution for the variance \code{sigma2} of the Gaussian process. \bold{Warning:} if a low-rank approximation is used, then \code{sigma2} corresponds to the variance of the iid zero-mean Gaussian variables. Default is \code{NULL}.
##' @param log.prior.phi a function corresponding to the log-density of the prior distribution for the scale parameter of the Matern correlation function; default is \code{NULL}.
##' @param log.prior.nugget optional: a function corresponding to the log-density of the prior distribution for the variance of the nugget effect; default is \code{NULL} with no nugget incorporated in the model; default is \code{NULL}.
##' @param uniform.sigma2 a vector of length two, corresponding to the lower and upper limit of the uniform prior on \code{sigma2}. Default is \code{NULL}.
##' @param uniform.phi a vector of length two, corresponding to the lower and upper limit of the uniform prior on \code{phi}. Default is \code{NULL}.
##' @param uniform.nugget a vector of length two, corresponding to the lower and upper limit of the uniform prior on \code{tau2}. Default is \code{NULL}.
##' @param log.normal.sigma2 a vector of length two, corresponding to the mean and standard deviation of the distribution on the log scale for the log-normal prior on \code{sigma2}. Default is \code{NULL}.
##' @param log.normal.phi a vector of length two, corresponding to the mean and standard deviation of the distribution on the log scale for the log-normal prior on \code{phi}. Default is \code{NULL}.
##' @param log.normal.nugget a vector of length two, corresponding to the mean and standard deviation of the distribution on the log scale for the log-normal prior on \code{tau2}. Default is \code{NULL}.
##' @return a list corresponding the prior distributions for each model parameter.
##' @seealso See "Priors definition" in the Details section of the \code{\link{binomial.logistic.Bayes}} function.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
control.prior <- function(beta.mean,beta.covar, log.prior.sigma2=NULL,
                                    log.prior.phi=NULL,log.prior.nugget=NULL,
                                    uniform.sigma2=NULL,log.normal.sigma2=NULL, 
                                    uniform.phi=NULL,log.normal.phi=NULL,
                                    uniform.nugget=NULL,log.normal.nugget=NULL) { 
   if(dim(as.matrix(beta.covar))[1]*
      dim(as.matrix(beta.covar))[2] != (length(beta.mean)^2)) {
      stop("incompatible values for beta.covar and beta.mean")
   }
   if(length(log.prior.sigma2) >0 && length(log.prior.sigma2(runif(1)))!=1) stop("invalid prior for sigma2")      
   if(length(log.prior.phi) >0 && length(log.prior.phi(runif(1)))!=1) stop("invalid prior for phi")  
   if(length(log.prior.nugget)>0) {
   	if(length(log.prior.nugget(runif(1)))!=1) stop("invalid prior for the nugget effect")  
   }         
    if(prod(t(beta.covar)==beta.covar)==0) stop("beta.covar is not symmetric.")
    if(prod(eigen(beta.covar)$values) <= 10e-07) stop("beta.covar is not positive definite.")
    if(length(uniform.sigma2) > 0 && any(uniform.sigma2 < 0)) stop("uniform.sigma2 must be positive.")	
    if(length(uniform.phi) > 0 && any(uniform.phi < 0)) stop("uniform.phi must be positive.")	
    if(length(uniform.nugget) > 0 && any(uniform.nugget < 0)) stop("uniform.nugget must be positive.")	

    if(length(log.normal.sigma2) > 0 && log.normal.sigma2[2] < 0) stop("the second element of log.normal.sigma2 must be positive.")	
    if(length(log.normal.phi) > 0 && log.normal.phi[2] < 0) stop("the second element of log.normal.phi must be positive.")	
    if(length(log.normal.nugget) > 0 && log.normal.nugget[2] < 0) stop("the second element of log.normal.nugget must be positive.")

    if((length(log.prior.sigma2) > 0 & (length(uniform.sigma2)>0 | length(log.normal.sigma2))) | 
       ((length(log.prior.sigma2) > 0 | length(uniform.sigma2)>0) & length(log.normal.sigma2)) |
       ((length(log.prior.sigma2) > 0 | length(log.normal.sigma2)>0) & length(uniform.sigma2))) stop("only one prior for sigma2 must be provided.")
    if((length(log.prior.phi) > 0 & (length(uniform.phi)>0 | length(log.normal.phi))) | 
       ((length(log.prior.phi) > 0 | length(uniform.phi)>0) & length(log.normal.phi)) |
       ((length(log.prior.phi) > 0 | length(log.normal.phi)>0) & length(uniform.phi))) stop("only one prior for phi must be provided.")       
    if((length(log.prior.nugget) > 0 & (length(uniform.nugget)>0 | length(log.normal.nugget))) | 
       ((length(log.prior.nugget) > 0 | length(uniform.nugget)>0) & length(log.normal.nugget)) |
       ((length(log.prior.nugget) > 0 | length(log.normal.nugget)>0) & length(uniform.nugget))) stop("only one prior for the variance of the nugget effect must be provided.")     
    if(length(uniform.sigma2) > 0 && length(uniform.sigma2) !=2) stop("wrong length of log.normal.sigma2")
    if(length(uniform.phi) > 0 && length(uniform.phi) !=2) stop("wrong length of log.normal.phi")
    if(length(uniform.nugget) > 0 && length(uniform.nugget) !=2) stop("wrong length of log.normal.nugget")
    if(length(log.normal.sigma2) > 0 && length(log.normal.sigma2) !=2) stop("wrong length of log.normal.sigma2")
    if(length(log.normal.phi) > 0 && length(log.normal.phi) !=2) stop("wrong length of log.normal.phi")
    if(length(log.normal.nugget) > 0 && length(log.normal.nugget) !=2) stop("wrong length of log.normal.nugget")
    
    if(length(uniform.sigma2) > 0 && (uniform.sigma2[2] < uniform.sigma2[1])) stop("wrong value for uniform.sigma2: the value for the upper limit is smaller than the lower limit.")
    if(length(uniform.phi) > 0 && (uniform.phi[2] < uniform.phi[1])) stop("wrong value for uniform.phi: the value for the upper limit is smaller than the lower limit.")
    if(length(uniform.nugget) > 0 && (uniform.nugget[2] < uniform.nugget[1])) stop("wrong value for uniform.nugget: the value for the upper limit is smaller than the lower limit.")
   
       if(length(uniform.sigma2) > 0 && any(uniform.sigma2<0)) stop("uniform.sigma2 must be positive.")
    if(length(uniform.phi) > 0 && any(uniform.phi<0)) stop("uniform.phi must be positive.")
    if(length(uniform.nugget) > 0 && any(uniform.nugget<0)) stop("uniform.nugget must be positive.")
    
   if(length(uniform.sigma2) > 0) {
   	   log.prior.sigma2 <- function(sigma2) dunif(sigma2,uniform.sigma2[1], uniform.sigma2[2],log=TRUE)
   } else if(length(log.normal.sigma2) > 0) {
   	   log.prior.sigma2 <- function(sigma2) dlnorm(sigma2,meanlog=log.normal.sigma2[1], 
   	                                                             log.normal.sigma2[2],log=TRUE)
   }
   if(length(uniform.phi) > 0) {
   	   log.prior.phi <- function(phi) dunif(phi,uniform.phi[1], uniform.phi[2],log=TRUE)
   } else if(length(log.normal.phi) > 0) {
   	   log.prior.phi <- function(phi) dlnorm(phi,meanlog=log.normal.phi[1], 
   	                                                            log.normal.phi[2],log=TRUE)
   }
   if(length(uniform.nugget) > 0) {
   	   log.prior.nugget <- function(tau2) dunif(tau2,uniform.nugget[1], uniform.nugget[2],log=TRUE)
   } else if(length(log.normal.nugget) > 0) {
   	   log.prior.nugget <- function(tau2) dlnorm(tau2,meanlog=log.normal.nugget[1], 
   	                                                             log.normal.nugget[2],log=TRUE)
   }
   if(length(beta.mean)!=dim(as.matrix(beta.covar))[1] |
      length(beta.mean)!=dim(as.matrix(beta.covar))[2]) {
   	 stop("invalid mean or covariance matrix for the prior of beta")
   }
   out <- list(m=beta.mean,V.inv=solve(beta.covar),
   log.prior.phi=log.prior.phi,
   log.prior.sigma2=log.prior.sigma2,
   log.prior.nugget=log.prior.nugget) 
   return(out)
}

##' @title Control settings for the MCMC algorithm used for Bayesian inference
##' @description This function defines the different tuning parameter that are used in the MCMC algorithm for Bayesian inference.
##' @param n.sim total number of simulations.
##' @param burnin initial number of samples to be discarded.
##' @param thin value used to retain only evey \code{thin}-th sampled value.
##' @param h.theta1 starting value of the tuning parameter of the proposal distribution for \eqn{\theta_{1} = \log(\sigma^2)/2}. See 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param h.theta2 starting value of the tuning parameter of the proposal distribution for \eqn{\theta_{2} = \log(\sigma^2/\phi^{2 \kappa})}. See 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param h.theta3 starting value of the tuning parameter of the proposal distribution for \eqn{\theta_{3} = \log(\tau^2)}. See 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param L.S.lim an atomic value or a vector of length 2 that is used to define the number of steps used at each iteration in the Hamiltonian Monte Carlo algorithm to update the spatial random effect; if a single value is provided than the number of steps is kept fixed, otherwise if a vector of length 2 is provided the number of steps is simulated at each iteration as \code{floor(runif(1,L.S.lim[1],L.S.lim[2]+1))}.
##' @param epsilon.S.lim an atomic value or a vector of length 2 that is used to define the stepsize used at each iteration in the Hamiltonian Monte Carlo algorithm to update the spatial random effect; if a single value is provided than the stepsize is kept fixed, otherwise if a vector of length 2 is provided the stepsize is simulated at each iteration as \code{runif(1,epsilon.S.lim[1],epsilon.S.lim[2])}.
##' @param start.beta starting value for the regression coefficients \code{beta}.
##' @param start.sigma2 starting value for \code{sigma2}.
##' @param start.phi starting value for \code{phi}.
##' @param start.S starting value for the spatial random effect.
##' @param start.nugget starting value for the variance of the nugget effect; default is \code{NULL} if the nugget effect is not present.
##' @param c1.h.theta1 value of \eqn{c_{1}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(sigma2)/2}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param c2.h.theta1 value of \eqn{c_{2}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(sigma2)/2}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param c1.h.theta2 value of \eqn{c_{1}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(sigma2.curr/(phi.curr^(2*kappa)))}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param c2.h.theta2 value of \eqn{c_{2}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(sigma2.curr/(phi.curr^(2*kappa)))}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param c1.h.theta3 value of \eqn{c_{1}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(tau2)}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param c2.h.theta3 value of \eqn{c_{2}} used to adaptively tune the variance of the Gaussian proposal for the transformed parameter \code{log(tau2)}; see 'Details' in \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param linear.model logical; if 	\code{linear.model=TRUE}, the control parameters are set for the geostatistical linear model. Default is \code{linear.model=FALSE}.
##' @param binary logical; if \code{binary=TRUE}, the control parameters are set the binary geostatistical model. Default is \code{binary=FALSE}.
##' @return an object of class "mcmc.Bayes.PrevMap".
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export

control.mcmc.Bayes <- function(n.sim,burnin,thin,
                               h.theta1,h.theta2,h.theta3=NULL,
                               L.S.lim=NULL, epsilon.S.lim=NULL,
                               start.beta,start.sigma2,
                               start.phi,start.S,start.nugget=NULL,
                               c1.h.theta1=0.01,c2.h.theta1=0.0001,
                               c1.h.theta2=0.01,c2.h.theta2=0.0001,
                               c1.h.theta3=0.01,c2.h.theta3=0.0001,
                               linear.model=FALSE,binary=FALSE) {
   if(!linear.model & !binary) epsilon.S.lim <- as.numeric(epsilon.S.lim)
   if(!linear.model & !binary && (any(length(epsilon.S.lim)==c(1,2))==FALSE)) stop("epsilon.S.lim must be: atomic; a vector of length 2.")  
   
   if(!linear.model & !binary && (any(length(L.S.lim)==c(1,2))==FALSE))  stop("L.S.lim must be: atomic; a vector of length 2; or NULL for the linear Gaussian model.")                     
   
   	if(n.sim < burnin) stop("n.sim cannot be smaller than burnin.")
	
	if((n.sim-burnin)%%thin!=0) stop("thin must be a divisor of (n.sim-burnin)")
	
	if(length(start.nugget) > 0 & length(h.theta3) == 0) stop("h.theta3 is missing.")	
	
	if(length(start.nugget) > 0 & length(c1.h.theta3) == 0) stop("c1.h.theta3 is missing.")	
	
	if(length(start.nugget) > 0 & length(c2.h.theta3) == 0) stop("c2.h.theta3 is missing.")	
	
	if(length(start.nugget) == 0 & length(h.theta3) > 0) stop("start.nugget is missing.")
	
	if(c1.h.theta1 < 0) stop("c1.h.theta1 must be positive.")
	
	if(c1.h.theta2 < 0) stop("c1.h.theta2 must be positive.")
	
	if(length(c1.h.theta3) > 0 && c1.h.theta3 < 0) stop("c1.h.theta3 must be positive.")
	
	if(c2.h.theta1 < 0 | c2.h.theta1 > 1) stop("c2.h.theta1 must be between 0 and 1.")
	
	if(c2.h.theta2 < 0 | c2.h.theta2 > 1) stop("c2.h.theta2 must be between 0 and 1.")
	
	if(length(c2.h.theta3) > 0 && (c2.h.theta3 < 0 | c2.h.theta3 > 1)) stop("c2.h.theta3 must be between 0 and 1.")	
	
	if(!linear.model & !binary && (length(epsilon.S.lim)==2 & epsilon.S.lim[2] < epsilon.S.lim[1])) stop("The first element of epsilon.S.lim must be smaller than the second.")	
    
    if(!linear.model & !binary && (length(L.S.lim)==2 & L.S.lim[2] < epsilon.S.lim[1])) stop("The first element of L.S.lim must be smaller than the second.")	
    
    if(linear.model) {
         out <- list(n.sim=n.sim,burnin=burnin,thin=thin,
                   h.theta1=h.theta1,h.theta2=h.theta2,
                   h.theta3=h.theta3,
                   start.beta=start.beta,start.sigma2=start.sigma2,
                   start.phi=start.phi,start.nugget=start.nugget,
                   c1.h.theta1=c1.h.theta1,c2.h.theta1=c2.h.theta1,
                   c1.h.theta2=c1.h.theta2,c2.h.theta2=c2.h.theta2,
                   c1.h.theta3=c1.h.theta3,c2.h.theta3=c2.h.theta3)	
   } else if (binary) {
   	    if(length(start.nugget) != 0) warning("inclusion of the nugget effect is not available for the binary model") 
   	    out <- list(n.sim=n.sim,burnin=burnin,thin=thin,
                   h.theta1=h.theta1,h.theta2=h.theta2,
                   start.beta=start.beta,start.sigma2=start.sigma2,
                   start.phi=start.phi,
                   c1.h.theta1=c1.h.theta1,c2.h.theta1=c2.h.theta1,
                   c1.h.theta2=c1.h.theta2,c2.h.theta2=c2.h.theta2)	
   } else {
         out <- list(n.sim=n.sim,burnin=burnin,thin=thin,
                   h.theta1=h.theta1,h.theta2=h.theta2,
                   h.theta3=h.theta3,
                   L.S.lim=L.S.lim, epsilon.S.lim=epsilon.S.lim,
                   start.beta=start.beta,start.sigma2=start.sigma2,
                   start.phi=start.phi,
                   start.S=start.S,start.nugget=start.nugget,
                   c1.h.theta1=c1.h.theta1,c2.h.theta1=c2.h.theta1,
                   c1.h.theta2=c1.h.theta2,c2.h.theta2=c2.h.theta2,
                   c1.h.theta3=c1.h.theta3,c2.h.theta3=c2.h.theta3)	
   }
   class(out) <- "mcmc.Bayes.PrevMap"                
   out
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom geoR varcov.spatial
binomial.geo.Bayes <- function(formula,units.m,coords,data,
                                                  ID.coords,
                                                  control.prior,
                                                  control.mcmc,
                                                  kappa,messages) {
   kappa <- as.numeric(kappa)
   if(length(control.prior$log.prior.nugget)!=length(control.mcmc$start.nugget)) {
      stop("missing prior or missing starting value for the nugget effect")	
   }
   if(any(is.na(data))) stop("missing values are not accepted")                                                                                             
   if(length(control.mcmc$L.S.lim)==0) stop("L.S.lim must be provided in control.mcmc")
   if(length(control.mcmc$epsilon.S.lim)==0) stop("epsilon.S.lim must be provided in control.mcmc")
   	coords <- as.matrix(model.frame(coords,data))
    units.m <-  as.numeric(model.frame(units.m,data)[,1])   
    if(is.numeric(units.m)==FALSE) stop("'units.m' must be numeric.")
	if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")
   out <- list()
   if(length(ID.coords)==0) {
   	   mf <- model.frame(formula,data=data)
	   y <- as.numeric(model.response(mf))
	   n <- length(y)
	   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))	   
   	   p <- ncol(D)
   	   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   	   U <- dist(coords)
   	   V.inv <- control.prior$V.inv
   	   m <- control.prior$m   	  
   	   nu <- 2*kappa
       log.prior.theta1_2 <- function(theta1,theta2) {
         sigma2 <- exp(2*theta1)
         phi <- (exp(2*theta1-theta2))^(1/nu)
         log(phi)+log(sigma2)+control.prior$log.prior.sigma2(sigma2)+
         control.prior$log.prior.phi(phi)	
       }
      
       log.prior.theta3 <- function(theta3) {
      	 tau2 <- exp(theta3)
      	 theta3+control.prior$log.prior.nugget(tau2)
       }
      
       grad.log.posterior.S <- function(S,mu,sigma2,param.post) {
	     h <- units.m*exp(S)/(1+exp(S))
	     as.numeric(-param.post$R.inv%*%(S-mu)/sigma2+y-h)
       } 
       
       fixed.nugget <- length(control.mcmc$start.nugget)==0
       
       log.posterior <- function(theta1,theta2,beta,S,param.post,theta3=NULL) {
	     sigma2 <- exp(2*theta1)
	     mu <- D%*%beta
	     diff.beta <- beta-m
	     diff.S <- S-mu
	     if(length(theta3)==1) {
	     	lp.theta3 <- log.prior.theta3(theta3)
	     } else {
	     	lp.theta3 <- 0
	     }	
	     out <- log.prior.theta1_2(theta1,theta2)+lp.theta3+
	     -0.5*(p*log(sigma2)+t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
	     -0.5*(n*log(sigma2)+param.post$ldetR+
	     t(diff.S)%*%param.post$R.inv%*%(diff.S)/sigma2)+
         sum(y*S-units.m*log(1+exp(S)))
         as.numeric(out)
      }
      
      acc.theta1 <- 0
      acc.theta2 <- 0
      if(fixed.nugget==FALSE) acc.theta3 <- 0
      acc.S <- 0 
      n.sim <- control.mcmc$n.sim 
      burnin <- control.mcmc$burnin
      thin <- control.mcmc$thin
      h.theta1 <- control.mcmc$h.theta1
      h.theta2 <- control.mcmc$h.theta2
      h.theta3 <- control.mcmc$h.theta3
      
      c1.h.theta1 <- control.mcmc$c1.h.theta1
      c2.h.theta1 <- control.mcmc$c2.h.theta1
      c1.h.theta2 <- control.mcmc$c1.h.theta2
      c2.h.theta2 <- control.mcmc$c2.h.theta2
      c1.h.theta3 <- control.mcmc$c1.h.theta3
      c2.h.theta3 <- control.mcmc$c2.h.theta3
      
      epsilon.S.lim <-  control.mcmc$epsilon.S.lim
      L.S.lim <- control.mcmc$L.S.lim
      if(!fixed.nugget) {
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+3)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi","tau^2")               
         tau2.curr <- control.mcmc$start.nugget
         theta3.curr <- log(tau2.curr)
      } else {  
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")    
         theta3.curr <- NULL
	  	 tau2.curr <- 0
      }   
      out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=n)

      
      # Initialise all the parameters
      sigma2.curr <- control.mcmc$start.sigma2
	  phi.curr <- control.mcmc$start.phi	  
	  theta1.curr <- 0.5*log(sigma2.curr)
	  theta2.curr <- log(sigma2.curr/(phi.curr^nu))
	  beta.curr <- control.mcmc$start.beta
	  S.curr <- control.mcmc$start.S
	
	  # Compute log-posterior density
	  R.curr <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                        cov.pars=c(1,phi.curr),
	                                        kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov

	  param.post.curr <- list()
	  param.post.curr$R.inv <- solve(R.curr)
	  param.post.curr$ldetR <- determinant(R.curr)$modulus
	  lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,
	                                    param.post.curr,theta3.curr)
      h1 <- rep(NA,n.sim) 
      h2 <- rep(NA,n.sim) 
      if(!fixed.nugget) h3 <- rep(NA,n.sim)              
      for(i in 1:n.sim) {
      	
         # Update theta1 
    	 theta1.prop <- theta1.curr+h.theta1*rnorm(1)
    	 sigma2.prop <- exp(2*theta1.prop)
	     phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
    	 R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.prop)$varcov 
         
         param.post.prop <- list()	     
	     param.post.prop$R.inv <- solve(R.prop)
	     param.post.prop$ldetR <- determinant(R.prop)$modulus
	     lp.prop <- log.posterior(theta1.prop,theta2.curr,beta.curr,S.curr,
	                                          param.post.prop,theta3.curr) 	                                         	     
	     if(log(runif(1)) < lp.prop-lp.curr) {
	        theta1.curr <- theta1.prop
	        param.post.curr <- param.post.prop
	        lp.curr <- lp.prop
	        sigma2.curr <- sigma2.prop
	        phi.curr <- phi.prop
	        acc.theta1 <- acc.theta1+1
	     } 
	     rm(theta1.prop,sigma2.prop,phi.prop)   	                                          
         h1[i] <- h.theta1 <- max(0,h.theta1 + 
	                       (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                   
    	     
    	 # Update theta2 
    	 theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	     phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
    	 R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov 
         
         param.post.prop <- list()	     
	     param.post.prop$R.inv <- solve(R.prop)
	     param.post.prop$ldetR <- determinant(R.prop)$modulus
	     lp.prop <- log.posterior(theta1.curr,theta2.prop,beta.curr,S.curr,
	                                          param.post.prop,theta3.curr) 	                                         
	     
	     if(log(runif(1)) < lp.prop-lp.curr) {
	        theta2.curr <- theta2.prop
	        param.post.curr <- param.post.prop
	        lp.curr <- lp.prop
	        phi.curr <- phi.prop
	        acc.theta2 <- acc.theta2+1
	     }    	                      
	     rm(theta2.prop,phi.prop)                     
         h2[i] <- h.theta2 <- max(0,h.theta2 + 
	                       (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45)) 
	                                            
         # Update theta3 
         if(!fixed.nugget) { 
    	    theta3.prop <- theta3.curr+h.theta3*rnorm(1)
	        tau2.prop <- exp(theta3.prop)
    	    R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.curr),
	                                         kappa=kappa,nugget=tau2.prop/sigma2.curr)$varcov 
         
            param.post.prop <- list()	     
	        param.post.prop$R.inv <- solve(R.prop)
	        param.post.prop$ldetR <- determinant(R.prop)$modulus
	        lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,
	                                          param.post.prop,theta3.prop) 	                                         
	     
	        if(log(runif(1)) < lp.prop-lp.curr) {
	           theta3.curr <- theta3.prop
	           param.post.curr <- param.post.prop
	           lp.curr <- lp.prop
	           tau2.curr <- tau2.prop
	           acc.theta3 <- acc.theta3+1
	        }    	                      
	        rm(theta3.prop,tau2.prop)                     
            h3[i] <- h.theta3 <- max(0,h.theta3 + 
	                       (c1.h.theta3*i^(-c2.h.theta3))*(acc.theta3/i-0.45)) 
	     }               
	                          
         # Update beta
         A <- t(D)%*%param.post.curr$R.inv
         cov.beta <- solve(V.inv+A%*%D)
         mean.beta <- as.numeric(cov.beta%*%(V.inv%*%m+A%*%S.curr))
         beta.curr <- as.numeric(mean.beta+sqrt(sigma2.curr)*
                             t(chol(cov.beta))%*%rnorm(p))
         lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,
	                                                 param.post.curr,theta3.curr)
         mu.curr <- D%*%beta.curr	   
      
         # Update S
         if(length(epsilon.S.lim)==2) {
            epsilon.S <- runif(1,epsilon.S.lim[1],epsilon.S.lim[2])	
         } else {
         	epsilon.S <- epsilon.S.lim
         }    
	     q.S <- S.curr
	     p.S = rnorm(n) 
         current_p.S = p.S
         p.S = p.S + epsilon.S * 
                     grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)/ 2

         if(length(L.S.lim)==2) {
            L.S <- floor(runif(1,L.S.lim[1],L.S.lim[2]+1))            
         } else {
         	L.S <- L.S.lim
         }                        
         for (j in 1:L.S) {
           q.S = q.S + epsilon.S * p.S
           if (j!=L.S) {p.S = p.S + epsilon.S *  
           grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)
           }
         }
         if(any(is.nan(q.S))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")

         p.S = p.S + epsilon.S* 
                  grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)/2
         p.S = -p.S
         current_U = -lp.curr
         current_K = sum(current_p.S^2) / 2
         proposed_U =  -log.posterior(theta1.curr,theta2.curr,beta.curr,q.S,
	                                    param.post.curr,theta3.curr)
         proposed_K = sum(p.S^2) / 2
         if(any(is.na(proposed_U) |  is.na(proposed_K))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")
         if (log(runif(1)) < current_U-proposed_U+current_K-proposed_K) {
            S.curr <- q.S  # accept
            lp.curr <- -proposed_U
         }
      
         if(i > burnin & (i-burnin)%%thin==0) {
		    out$S[(i-burnin)/thin,] <- S.curr	
		    if(fixed.nugget) {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr)
		    	} else {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr,tau2.curr)
		    }
   	     }
         if(messages) cat("Iteration",i,"out of",n.sim,"\r")
         flush.console()                            
      }  	
   } else if(length(ID.coords)>0) {
   	   coords <- as.matrix(unique(coords))
   	   mf <- model.frame(formula,data=data)
	   y <- as.numeric(model.response(mf))
	   n <- length(y)
	   n.x <- nrow(coords)
	   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
   	   p <- ncol(D)
   	   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   	   U <- dist(coords)
   	   V.inv <- control.prior$V.inv
   	   m <- control.prior$m   	  
   	   nu <- 2*kappa
   	   C.S <- t(sapply(1:n.x,function(i) ID.coords==i))
       log.prior.theta1_2 <- function(theta1,theta2) {
         sigma2 <- exp(2*theta1)
         phi <- (exp(2*theta1-theta2))^(1/nu)
         log(phi)+log(sigma2)+control.prior$log.prior.sigma2(sigma2)+
         control.prior$log.prior.phi(phi)	
       }
      
       log.prior.theta3 <- function(theta3) {
      	 tau2 <- exp(theta3)
      	 theta3+control.prior$log.prior.nugget(tau2)
       }
       
       fixed.nugget <- length(control.mcmc$start.nugget)==0
       
       log.posterior <- function(theta1,theta2,beta,S,param.post,theta3=NULL) {
	     sigma2 <- exp(2*theta1)
	     mu <- as.numeric(D%*%beta)
	     diff.beta <- beta-m
	     eta <- mu+S[ID.coords]
	     if(length(theta3)==1) {
	     	lp.theta3 <- log.prior.theta3(theta3)
	     } else {
	     	lp.theta3 <- 0
	     }	
	     as.numeric(log.prior.theta1_2(theta1,theta2)+lp.theta3+
	     -0.5*(p*log(sigma2)+t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
	     -0.5*(n.x*log(sigma2)+param.post$ldetR+
	     t(S)%*%param.post$R.inv%*%(S)/sigma2)+
         sum(y*eta-units.m*log(1+exp(eta))))
      }
      
      integrand <- function(beta,sigma2,S) {
      	 diff.beta <- beta-m
      	 eta <- as.numeric(D%*%beta+S[ID.coords])
      	 -0.5*(t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
      	 sum(y*eta-units.m*log(1+exp(eta)))
      }
      
      grad.integrand <- function(beta,sigma2,S) {
      	 eta <- as.numeric(D%*%beta+S[ID.coords])
      	 h <- units.m*exp(eta)/(1+exp(eta))
      	 as.numeric(-V.inv%*%(beta-m)/sigma2+t(D)%*%(y-h))
      }
      
      hessian.integrand <- function(beta,sigma2,S) {
      	 eta <- as.numeric(D%*%beta+S[ID.coords])
      	 h1 <- units.m*exp(eta)/((1+exp(eta))^2)
      	 -V.inv/sigma2-t(D)%*%(D*h1)
      }
      
      grad.log.posterior.S <- function(S,mu,sigma2,param.post) {
      	 eta <- mu+S[ID.coords]     	
	     h <- units.m*exp(eta)/(1+exp(eta))
	     as.numeric(-param.post$R.inv%*%S/sigma2+ 
	     sapply(1:n.x,function(i) sum((y-h)[C.S[i,]])))
      } 
      
      acc.theta1 <- 0
      acc.theta2 <- 0
      if(fixed.nugget==FALSE) acc.theta3 <- 0
      acc.S <- 0 
      n.sim <- control.mcmc$n.sim 
      burnin <- control.mcmc$burnin
      thin <- control.mcmc$thin
      h.theta1 <- control.mcmc$h.theta1
      h.theta2 <- control.mcmc$h.theta2
      h.theta3 <- control.mcmc$h.theta3
      
      c1.h.theta1 <- control.mcmc$c1.h.theta1
      c2.h.theta1 <- control.mcmc$c2.h.theta1
      c1.h.theta2 <- control.mcmc$c1.h.theta2
      c2.h.theta2 <- control.mcmc$c2.h.theta2
      c1.h.theta3 <- control.mcmc$c1.h.theta3
      c2.h.theta3 <- control.mcmc$c2.h.theta3
      
      L.S.lim <- control.mcmc$L.S.lim
      epsilon.S.lim <-  control.mcmc$epsilon.S.lim
      out <- list()
      if(!fixed.nugget) {
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+3)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi","tau^2")               
         tau2.curr <- control.mcmc$start.nugget
         theta3.curr <- log(tau2.curr)
      } else {  
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")    
         theta3.curr <- NULL
	  	 tau2.curr <- 0
      }   
      out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=n.x)

      # Initialise all the parameters
      sigma2.curr <- control.mcmc$start.sigma2
	  phi.curr <- control.mcmc$start.phi
	  tau2.curr <- control.mcmc$start.nugget
	  theta1.curr <- 0.5*log(sigma2.curr)
	  theta2.curr <- log(sigma2.curr/(phi.curr^nu))
	  if(fixed.nugget==FALSE) {
	  	theta3.curr <- log(tau2.curr)
	  } else {
	  	theta3.curr <- NULL
	  	tau2.curr <- 0
	  }	
	  beta.curr <- control.mcmc$start.beta
	  S.curr <- control.mcmc$start.S
	  
	  # Compute the log-posterior density
	  if(fixed.nugget) {
	      R.curr <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                            cov.pars=c(1,phi.curr),
	                                            kappa=kappa)$varcov	   	
	  } else {
	   	   R.curr <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                            cov.pars=c(1,phi.curr),
	                                            kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov
	  }

	  param.post.curr <- list()
	  param.post.curr$R.inv <- solve(R.curr)
	  param.post.curr$ldetR <- determinant(R.curr)$modulus
	  lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,
	                                    param.post.curr,theta3.curr)
      mu.curr <- D%*%beta.curr
       
      h1 <- rep(NA,n.sim) 
      h2 <- rep(NA,n.sim) 
      if(!fixed.nugget) h3 <- rep(NA,n.sim)   	                
      for(i in 1:n.sim) {
      	
         # Update theta1 
    	 theta1.prop <- theta1.curr+h.theta1*rnorm(1)
    	 sigma2.prop <- exp(2*theta1.prop)
	     phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
    	 R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.prop)$varcov 
       
         param.post.prop <- list()	     
	     param.post.prop$R.inv <- solve(R.prop)
	     param.post.prop$ldetR <- determinant(R.prop)$modulus
	     lp.prop <- log.posterior(theta1.prop,theta2.curr,beta.curr,S.curr,
	                                          param.post.prop,theta3.curr) 	                                         	     
	     if(log(runif(1)) < lp.prop-lp.curr) {
	        theta1.curr <- theta1.prop
	        param.post.curr <- param.post.prop
	        lp.curr <- lp.prop
	        sigma2.curr <- sigma2.prop
	        phi.curr <- phi.prop
	        acc.theta1 <- acc.theta1+1
	     } 
	     rm(theta1.prop,sigma2.prop,phi.prop)   	                                          
         h1[i] <- h.theta1 <- max(0,h.theta1 + 
	                       (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                  
    	     
    	  # Update theta2 
    	  theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	      phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
    	  R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov 
         
          param.post.prop <- list()	     
	      param.post.prop$R.inv <- solve(R.prop)
	      param.post.prop$ldetR <- determinant(R.prop)$modulus
	      lp.prop <- log.posterior(theta1.curr,theta2.prop,beta.curr,S.curr,
	                                           param.post.prop,theta3.curr) 	                                         
	     
	      if(log(runif(1)) < lp.prop-lp.curr) {
	         theta2.curr <- theta2.prop
	         param.post.curr <- param.post.prop
	         lp.curr <- lp.prop
	         phi.curr <- phi.prop
	         acc.theta2 <- acc.theta2+1
	      }    	                      
	      rm(theta2.prop,phi.prop)                     
          h2[i] <- h.theta2 <- max(0,h.theta2 + 
	                        (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))  
	                                            
          # Update theta3
          if(!fixed.nugget) { 
    	      theta3.prop <- theta3.curr+h.theta3*rnorm(1)
	          tau2.prop <- exp(theta3.prop)
    	      R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.curr),
	                                         kappa=kappa,nugget=tau2.prop/sigma2.curr)$varcov 
         
              param.post.prop <- list()	     
	          param.post.prop$R.inv <- solve(R.prop)
	          param.post.prop$ldetR <- determinant(R.prop)$modulus
	          lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,
	                                          param.post.prop,theta3.prop) 	                                         
	     
	          if(log(runif(1)) < lp.prop-lp.curr) {
	             theta3.curr <- theta3.prop
	             param.post.curr <- param.post.prop
	             lp.curr <- lp.prop
	             tau2.curr <- tau2.prop
	             acc.theta3 <- acc.theta3+1
	          }    	                      
	          rm(theta3.prop,tau2.prop)                     
              h3[i] <- h.theta3 <- max(0,h.theta3 + 
	                       (c1.h.theta3*i^(-c2.h.theta3))*(acc.theta3/i-0.45)) 
	     }               
	                          
         # Update beta
         max.beta <- maxBFGS(function(x) integrand(x,sigma2.curr,S.curr),
                                               function(x) grad.integrand(x,sigma2.curr,S.curr),
                                               function(x)  hessian.integrand(x,sigma2.curr,S.curr),
                                               start=beta.curr)
         Sigma.tilde <- solve(-max.beta$hessian)
         Sigma.tilde.inv <- -max.beta$hessian
         beta.prop <- as.numeric(max.beta$estimate+t(chol(Sigma.tilde))%*%rt(p,10))
         diff.b.prop <- beta.curr-max.beta$estimate
         diff.b.curr <- beta.prop-max.beta$estimate
         q.prop <- as.numeric(-0.5*t(diff.b.prop)%*%Sigma.tilde.inv%*%(diff.b.prop))
         q.curr <- as.numeric(-0.5*t(diff.b.curr)%*%Sigma.tilde.inv%*%(diff.b.curr))
       
         lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.prop,S.curr,
	                                          param.post.curr,theta3.curr) 
         if(log(runif(1)) < lp.prop+q.prop-lp.curr-q.curr) {
         	lp.curr <- lp.prop
         	beta.curr <- beta.prop
         	mu.curr <- D%*%beta.curr	   
         } 
                         
         # Update S
         if(length(epsilon.S.lim)==2) {
            epsilon.S <- runif(1,epsilon.S.lim[1],epsilon.S.lim[2])	
         } else {
         	epsilon.S <- epsilon.S.lim
         }    
         
         if(length(L.S.lim)==2) {
            L.S <- floor(runif(1,L.S.lim[1],L.S.lim[2]+1))            
         } else {
         	L.S <- L.S.lim
         }    
         
	     q.S <- S.curr
	     p.S = rnorm(n.x) 
         current_p.S = p.S
         p.S = p.S + epsilon.S * 
                     grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)/ 2
         for (j in 1:L.S) {
               q.S = q.S + epsilon.S * p.S
               if (j!=L.S) {p.S = p.S + epsilon.S *  
               grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)
            }
         }
         if(any(is.nan(q.S))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")

         p.S = p.S + epsilon.S* 
                  grad.log.posterior.S(q.S,mu.curr,sigma2.curr,param.post.curr)/2
         p.S = -p.S
         current_U = -lp.curr
         current_K = sum(current_p.S^2) / 2
         proposed_U =  -log.posterior(theta1.curr,theta2.curr,beta.curr,q.S,
	                                    param.post.curr,theta3.curr)
         proposed_K = sum(p.S^2) / 2
         if(any(is.na(proposed_U) |  is.na(proposed_K))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")
         if (log(runif(1)) < current_U-proposed_U+current_K-proposed_K) {
            S.curr <- q.S  # accept
            lp.curr <- -proposed_U
         }
      
         if(i > burnin & (i-burnin)%%thin==0) {
		    out$S[(i-burnin)/thin,] <- S.curr
		    if(fixed.nugget) {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr)
		    	} else {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr,tau2.curr)
		    	}		
   	        }
         if(messages) cat("Iteration",i,"out of",n.sim,"\r")
         flush.console()                            
      }  	
   }
   class(out) <- "Bayes.PrevMap"	
   out$y <- y
   out$units.m <- units.m
   out$D <- D
   out$coords <- coords
   out$kappa <- kappa
   out$ID.coords <- ID.coords  
   out$h1 <- h1
   out$h2 <- h2
   if(!fixed.nugget) out$h3 <- h3
   return(out)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
binomial.geo.Bayes.lr <- function(formula,units.m,coords,data,knots,
                                                       control.mcmc,control.prior,kappa,
                                                       messages) {	
       knots <- as.matrix(knots)
       coords <- as.matrix(model.frame(coords,data))
       units.m <-  as.numeric(model.frame(units.m,data)[,1])   
       if(is.numeric(units.m)==FALSE) stop("'units.m' must be numeric.")
	   if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")                                                       	                                                       	
   	   mf <- model.frame(formula,data=data)
	   y <- as.numeric(model.response(mf))
	   n <- length(y)
	   N <- nrow(knots)
	   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))	   
   	   p <- ncol(D)
   	   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   	   U.k <- as.matrix(pdist(coords,knots))
   	   V.inv <- control.prior$V.inv
   	   m <- control.prior$m   	  
   	   nu <- 2*kappa

       log.prior.theta1_2 <- function(theta1,theta2) {
         sigma2 <- exp(2*theta1)
         phi <- (exp(2*theta1-theta2))^(1/nu)
         log(phi/kappa)+control.prior$log.prior.sigma2(sigma2)+
         control.prior$log.prior.phi(phi)	
       }
             
        log.posterior <- function(theta1,theta2,beta,S,W) {
	     sigma2 <- exp(2*theta1)
	     mu <- as.numeric(D%*%beta)
	     diff.beta <- beta-m
	     eta <- as.numeric(mu+W)
	     as.numeric(log.prior.theta1_2(theta1,theta2)+
	     -0.5*(p*log(sigma2)+t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
	     -0.5*(N*log(sigma2)+sum(S^2)/sigma2)+
         sum(y*eta-units.m*log(1+exp(eta))))
       }
      
       integrand <- function(beta,sigma2,W) {
      	 diff.beta <- beta-m
      	 eta <- as.numeric(D%*%beta+W)
      	 -0.5*(t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
      	 sum(y*eta-units.m*log(1+exp(eta)))
       }
      
       grad.integrand <- function(beta,sigma2,W) {
       	 eta <- as.numeric(D%*%beta+W)
       	 h <- units.m*exp(eta)/(1+exp(eta))
       	 as.numeric(-V.inv%*%(beta-m)/sigma2+t(D)%*%(y-h))
       } 
      
       hessian.integrand <- function(beta,sigma2,W) {
      	 eta <- as.numeric(D%*%beta+W)
       	 h1 <- units.m*exp(eta)/((1+exp(eta))^2)
       	 -V.inv/sigma2-t(D)%*%(D*h1)
       }
      
       grad.log.posterior.S <- function(S,mu,sigma2,K) {
        	  eta <- as.numeric(mu+K%*%S)
	      h <- units.m*exp(eta)/(1+exp(eta))
	      as.numeric(-S/sigma2+t(K)%*%(y-h))
       } 
      
       c1.h.theta1 <- control.mcmc$c1.h.theta1
       c1.h.theta2 <- control.mcmc$c1.h.theta2
       
       c2.h.theta1 <- control.mcmc$c2.h.theta1
       c2.h.theta2 <- control.mcmc$c2.h.theta2
       
       h.theta1 <- control.mcmc$h.theta1
       h.theta2 <- control.mcmc$h.theta2
       
       acc.theta1 <- 0
       acc.theta2 <- 0
       acc.S <- 0 
       n.sim <- control.mcmc$n.sim 
       burnin <- control.mcmc$burnin
       thin <- control.mcmc$thin
       L.S.lim <- control.mcmc$L.S.lim
       epsilon.S.lim <-  control.mcmc$epsilon.S.lim
       out <- list()
       out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
       colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")
       theta3.curr <- NULL
	   tau2.curr <- 0
	  	 
       out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=N)
       out$const.sigma2 <- rep(NA,(n.sim-burnin)/thin)
       # Initialise all the parameters
       sigma2.curr <- control.mcmc$start.sigma2
	   phi.curr <- control.mcmc$start.phi
	   rho.curr <- phi.curr*2*sqrt(kappa)
	   theta1.curr <- 0.5*log(sigma2.curr)
	   theta2.curr <- log(sigma2.curr/(phi.curr^nu))

	   beta.curr <- control.mcmc$start.beta
	   S.curr <- control.mcmc$start.S
	
	   # Compute the log-posterior density
	   K.curr <- matern.kernel(U.k,rho.curr,kappa)	 
       W.curr <- as.numeric(K.curr%*%S.curr)

	   lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,S.curr,W.curr)
       mu.curr <- as.numeric(D%*%beta.curr	)
       
       h1 <- rep(NA,n.sim)
       h2 <- rep(NA,n.sim)
       for(i in 1:n.sim) {
      	
         # Update theta1 
    	 theta1.prop <- theta1.curr+h.theta1*rnorm(1)
    	 sigma2.prop <- exp(2*theta1.prop)
	     phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
	     rho.prop <- phi.prop*2*sqrt(kappa)
    	 K.prop <- matern.kernel(U.k,rho.curr,kappa)
         
         W.prop <- as.numeric(K.prop%*%S.curr)
	     lp.prop <- log.posterior(theta1.prop,theta2.curr,beta.curr,S.curr,W.prop) 	                                         	     
	     if(log(runif(1)) < lp.prop-lp.curr) {
	        theta1.curr <- theta1.prop
	        lp.curr <- lp.prop
	        sigma2.curr <- sigma2.prop
	        phi.curr <- phi.prop
	        rho.curr <- rho.prop
	        K.curr <- K.prop
	        W.curr <- W.prop
	        acc.theta1 <- acc.theta1+1
	     } 
	     rm(theta1.prop,sigma2.prop,phi.prop,K.prop,W.prop,rho.prop)   	                                          
         h1[i] <- h.theta1 <- max(0,h.theta1 + 
	                       (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                    
    	     
    	  # Update theta2 
    	  theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	      phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
	      rho.prop <- 2*phi.prop*sqrt(kappa)
    	  K.prop <- matern.kernel(U.k,rho.prop,kappa) 
         
         W.prop <- as.numeric(K.prop%*%S.curr)
	     lp.prop <- log.posterior(theta1.curr,theta2.prop,beta.curr,S.curr,W.prop) 	                                         
	     
	     if(log(runif(1)) < lp.prop-lp.curr) {
	        theta2.curr <- theta2.prop
	        lp.curr <- lp.prop
	        phi.curr <- phi.prop
	        rho.curr <- rho.prop
	        K.curr <- K.prop
	        W.curr <- W.prop
	        acc.theta2 <- acc.theta2+1
	     }    	                      
	     rm(theta2.prop,phi.prop,rho.prop,K.prop,W.prop)                     
         h2[i] <- h.theta2 <- max(0,h.theta2 + 
	                       (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))                                                          
	                          
         # Update beta
         max.beta <- maxBFGS(function(x) integrand(x,sigma2.curr,W.curr),
                                               function(x) grad.integrand(x,sigma2.curr,W.curr),
                                               function(x) hessian.integrand(x,sigma2.curr,W.curr),
                                               start=beta.curr)
         Sigma.tilde <- solve(-max.beta$hessian)
         Sigma.tilde.inv <- -max.beta$hessian
         beta.prop <- as.numeric(max.beta$estimate+t(chol(Sigma.tilde))%*%rt(p,10))
         diff.b.prop <- beta.curr-max.beta$estimate
         diff.b.curr <- beta.prop-max.beta$estimate
         q.prop <- as.numeric(-0.5*t(diff.b.prop)%*%Sigma.tilde.inv%*%(diff.b.prop))
         q.curr <- as.numeric(-0.5*t(diff.b.curr)%*%Sigma.tilde.inv%*%(diff.b.curr))
       
         lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.prop,S.curr,W.curr) 
         if(log(runif(1)) < lp.prop+q.prop-lp.curr-q.curr) {
         	lp.curr <- lp.prop
         	beta.curr <- beta.prop
         	mu.curr <- D%*%beta.curr	   
         } 
                        
         # Update S
         if(length(epsilon.S.lim)==2) {
            epsilon.S <- runif(1,epsilon.S.lim[1],epsilon.S.lim[2])	
         } else {
         	epsilon.S <- epsilon.S.lim
         }    
         
         if(length(L.S.lim)==2) {
            L.S <- floor(runif(1,L.S.lim[1],L.S.lim[2]+1))            
          } else {
         	L.S <- L.S.lim
          }    
	     q.S <- S.curr
	     p.S = rnorm(N) 
         current_p.S = p.S
         p.S = p.S + epsilon.S * 
                     grad.log.posterior.S(q.S,mu.curr,sigma2.curr,K.curr)/ 2
         for (j in 1:L.S) {
             q.S = q.S + epsilon.S * p.S
             if (j!=L.S) {p.S = p.S + epsilon.S *  
               grad.log.posterior.S(q.S,mu.curr,sigma2.curr,K.curr)
             }
         }
         if(any(is.nan(q.S))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")

         p.S = p.S + epsilon.S* 
                  grad.log.posterior.S(q.S,mu.curr,sigma2.curr,K.curr)/2
         p.S = -p.S
         current_U = -lp.curr
         current_K = sum(current_p.S^2) / 2
         W.prop <- K.curr%*%q.S
         proposed_U =  -log.posterior(theta1.curr,theta2.curr,beta.curr,q.S,W.prop)
         proposed_K = sum(p.S^2) / 2
         if(any(is.na(proposed_U) |  is.na(proposed_K))) stop("extremely large value (in absolute value) proposed for the spatial random effect; try the following actions: \n
         1) smaller values for epsilon.S.lim and/or L.S.lim; \n
         2) less diffuse prior for the covariance parameters.")
         if (log(runif(1)) < current_U-proposed_U+current_K-proposed_K) {
            S.curr <- q.S  # accept
            W.curr <- W.prop
            lp.curr <- -proposed_U
         }
      
         if(i > burnin & (i-burnin)%%thin==0) {
		    out$S[(i-burnin)/thin,] <- S.curr
		    out$estimate[(i-burnin)/thin,]  <- c(beta.curr,sigma2.curr,phi.curr)	
		    out$const.sigma2[(i-burnin)/thin] <- mean(apply(
		                                                   K.curr,1,function(r) sqrt(sum(r^2))))
   	     }
         if(messages) cat("Iteration",i,"out of",n.sim,"\r")
         flush.console()                            
      }
      class(out) <- "Bayes.PrevMap"	
      out$y <- y
      out$units.m <- units.m
      out$D <- D
      out$coords <- coords       
      out$kappa <- kappa  
      out$knots <- knots   
      out$h1 <- h1
      out$h2 <- h2
      return(out)
}

##' @title Bayesian estimation for the binomial logistic model
##' @description This function performs Bayesian estimation for a geostatistical binomial logistic model.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param units.m an object of class \code{\link{formula}} indicating the binomial denominators.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param ID.coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID.coords}}. These must be provided if, for example, spatial random effects are defined at household level but some of the covariates are at individual level. \bold{Warning}: the household coordinates must all be distinct otherwise see \code{\link{jitterDupCoords}}. Default is \code{NULL}.
##' @param control.prior output from \code{\link{control.prior}}.
##' @param control.mcmc output from \code{\link{control.mcmc.Bayes}}.
##' @param kappa value for the shape parameter of the Matern covariance function.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation is required. Default is \code{low.rank=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @details
##' This function performs Bayesian estimation for the parameters of the geostatistical binomial logistic model. Conditionally on a zero-mean stationary Gaussian process \eqn{S(x)} and mutually independent zero-mean Gaussian variables \eqn{Z} with variance \code{tau2}, the linear predictor assumes the form
##' \deqn{\log(p/(1-p)) = d'\beta + S(x) + Z,}
##' where \eqn{d} is a vector of covariates with associated regression coefficients \eqn{\beta}. The Gaussian process \eqn{S(x)} has isotropic Matern covariance function (see \code{\link{matern}}) with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. 
##' 
##' \bold{Priors definition.} Priors can be defined through the function \code{\link{control.prior}}. The hierarchical structure of the priors is the following. Let \eqn{\theta} be the vector of the covariance parameters \code{c(sigma2,phi,tau2)}; then each component of \eqn{\theta} has independent priors freely defined by the user. However, in  \code{\link{control.prior}} uniform and log-normal priors are also available as default priors for each of the covariance parameters. To remove the nugget effect \eqn{Z}, no prior should be defined for \code{tau2}. Conditionally on \code{sigma2}, the vector of regression coefficients \code{beta} has a multivariate Gaussian prior with mean \code{beta.mean} and covariance matrix \code{sigma2*beta.covar}, while in the low-rank approximation the covariance matrix is simply \code{beta.covar}. 
##' 
##' \bold{Updating the covariance parameters with a Metropolis-Hastings algorithm.} In the MCMC algorithm implemented in \code{binomial.logistic.Bayes}, the transformed parameters \deqn{(\theta_{1}, \theta_{2}, \theta_{3})=(\log(\sigma^2)/2,\log(\sigma^2/\phi^{2 \kappa}), \log(\tau^2))} are independently updated using a Metropolis Hastings algorithm. At the \eqn{i}-th iteration, a new value is proposed for each from a univariate Gaussian distrubion with variance \eqn{h_{i}^2} that is tuned using the following adaptive scheme \deqn{h_{i} = h_{i-1}+c_{1}i^{-c_{2}}(\alpha_{i}-0.45),} where \eqn{\alpha_{i}} is the acceptance rate at the \eqn{i}-th iteration, 0.45 is the optimal acceptance rate for a univariate Gaussian distribution, whilst \eqn{c_{1} > 0} and \eqn{0 < c_{2} < 1} are pre-defined constants. The starting values \eqn{h_{1}} for each of the parameters \eqn{\theta_{1}}, \eqn{\theta_{2}} and \eqn{\theta_{3}} can be set using the function \code{\link{control.mcmc.Bayes}} through the arguments \code{h.theta1}, \code{h.theta2} and \code{h.theta3}. To define values for \eqn{c_{1}} and \eqn{c_{2}}, see the documentation of \code{\link{control.mcmc.Bayes}}.
##' 
##' \bold{Hamiltonian Monte Carlo.} The MCMC algorithm in \code{binomial.logistic.Bayes} uses a Hamiltonian Monte Carlo (HMC) procedure to update the random effect \eqn{T=d'\beta + S(x) + Z}; see Neal (2011) for an introduction to HMC. HMC makes use of a postion vector, say \eqn{t}, representing the random effect \eqn{T}, and a momentum vector, say \eqn{q}, of the same length of the position vector, say \eqn{n}. Hamiltonian dynamics also have a physical interpretation where the states of the system are described by the position of a puck and its momentum (its mass times its velocity). The Hamiltonian function is then defined as a function of \eqn{t} and \eqn{q}, having the form \eqn{H(t,q) = -\log\{f(t | y, \beta, \theta)\} + q'q/2}, where \eqn{f(t | y, \beta, \theta)} is the conditional distribution of \eqn{T} given the data \eqn{y}, the regression parameters \eqn{\beta} and covariance parameters \eqn{\theta}. The system of Hamiltonian equations then defines the evolution of the system in time, which can be used to implement an algorithm for simulation from the posterior distribution of \eqn{T}. In order to implement the Hamiltonian dynamic on a computer, the Hamiltonian equations must be discretised. The \emph{leapfrog method} is then used for this purpose, where two tuning parameters should be defined: the stepsize \eqn{\epsilon} and the number of steps \eqn{L}. These respectively correspond to \code{epsilon.S.lim} and \code{L.S.lim} in the \code{\link{control.mcmc.Bayes}} function. However, it is advisable to let \eqn{epsilon} and \eqn{L} take different random values at each iteration of the HCM algorithm so as to account for the different variances amongst the components of the posterior of \eqn{T}. This can be done in \code{\link{control.mcmc.Bayes}} by defning \code{epsilon.S.lim} and \code{L.S.lim} as vectors of two elements, each of which represents the lower and upper limit of a uniform distribution used to generate values for  \code{epsilon.S.lim} and \code{L.S.lim}, respectively.
##' 
##' \bold{Using a two-level model to include household-level and individual-level information.}
##' When analysing data from household sruveys, some of the avilable information information might be at household-level (e.g. material of house, temperature) and some at individual-level (e.g. age, gender). In this case, the Gaussian spatial process \eqn{S(x)} and the nugget effect \eqn{Z} are defined at hosuehold-level in order to account for extra-binomial variation between and within households, respectively. 
##'
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), \code{sigma2} may take very different values from the actual variance of the Gaussian process to approximate. The function \code{\link{adjust.sigma2}} can then be used to (approximately) explore the range for \code{sigma2}. For example if the variance of the Gaussian process is \code{0.5}, then an approximate value for \code{sigma2} is \code{0.5/const.sigma2}, where \code{const.sigma2} is the value obtained with \code{\link{adjust.sigma2}}.
##' @return An object of class "Bayes.PrevMap".
##' The function \code{\link{summary.Bayes.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: matrix of the posterior samples of the model parameters.
##' @return \code{S}: matrix of the posterior samples for each component of the random effect. 
##' @return \code{const.sigma2}: vector of the values of the multiplicative factor used to adjust the values of \code{sigma2} in the low-rank approximation.
##' @return \code{y}: binomial observations.
##' @return \code{units.m}: binomial denominators.
##' @return \code{D}: matrix of covariarates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{kappa}: shape parameter of the Matern function.
##' @return \code{ID.coords}: set of ID values defined through the argument \code{ID.coords}.
##' @return \code{knots}: matrix of spatial knots used in the low-rank approximation.
##' @return \code{h1}: vector of values taken by the tuning parameter \code{h.theta1} at each iteration.
##' @return \code{h2}: vector of values taken by the tuning parameter \code{h.theta2} at each iteration.
##' @return \code{h3}: vector of values taken by the tuning parameter \code{h.theta3} at each iteration.
##' @return \code{call}: the matched call.
##' @seealso  \code{\link{control.mcmc.Bayes}},  \code{\link{control.prior}},\code{\link{summary.Bayes.PrevMap}}, \code{\link{matern}}, \code{\link{matern.kernel}}, \code{\link{create.ID.coords}}.
##' @references Neal, R. M. (2011) \emph{MCMC using Hamiltonian Dynamics}, In: Handbook of Markov Chain Monte Carlo (Chapter 5), Edited by Steve Brooks, Andrew Gelman, Galin Jones, and Xiao-Li Meng Chapman & Hall / CRC Press.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @examples
##' set.seed(1234)
##' data(data_sim)
##' # Select a subset of data_sim with 50 observations
##' n.subset <- 50
##' data_subset <- data_sim[sample(1:nrow(data_sim),n.subset),]
##' # Set the MCMC control parameters
##' control.mcmc <- control.mcmc.Bayes(n.sim=10,burnin=0,thin=1,
##'                            h.theta1=0.05,h.theta2=0.05,
##'                            L.S.lim=c(1,50),epsilon.S.lim=c(0.01,0.02),
##'                            start.beta=0,start.sigma2=1,start.phi=0.15,
##'                            start.S=rep(0,n.subset))
##' 
##' cp <- control.prior(beta.mean=0,beta.covar=1,
##'                              log.normal.phi=c(log(0.15),0.05),
##'                              log.normal.sigma2=c(log(1),0.1)) 
##' 
##' fit.Bayes <- binomial.logistic.Bayes(formula=y~1,coords=~x1+x2,units.m=~units.m,
##'                                      data=data_subset,control.prior=cp,
##'                                      control.mcmc=control.mcmc,kappa=2)  
##' summary(fit.Bayes)
##'                                     
##' par(mfrow=c(2,4))
##' autocor.plot(fit.Bayes,param="S",component.S="all")
##' autocor.plot(fit.Bayes,param="beta",component.beta=1)
##' autocor.plot(fit.Bayes,param="sigma2")
##' autocor.plot(fit.Bayes,param="phi")
##' trace.plot(fit.Bayes,param="S",component.S=30)
##' trace.plot(fit.Bayes,param="beta",component.beta=1)
##' trace.plot(fit.Bayes,param="sigma2")
##' trace.plot(fit.Bayes,param="phi")
##' @export
binomial.logistic.Bayes <- function(formula,units.m,coords,data,ID.coords=NULL,
                                                   control.prior,control.mcmc,kappa,low.rank=FALSE,
                                                   knots=NULL,messages=TRUE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.") 
    if(low.rank & length(ID.coords) > 0) stop("the low-rank approximation is not available for a two-levels model with logistic link function; see instead ?binary.probit.Bayes")
    if(class(control.mcmc)!="mcmc.Bayes.PrevMap") stop("control.mcmc must be of class 'mcmc.Bayes.PrevMap'")
    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(class(units.m)!="formula") stop("units.m must be a 'formula' object indicating the binomial denominators.")
    if(kappa < 0) stop("kappa must be positive.")
	if(!low.rank) {
		res <- binomial.geo.Bayes(formula=formula,units.m=units.m,coords=coords,
		           data=data,ID.coords=ID.coords,control.prior=control.prior,
		           control.mcmc=control.mcmc,
		           kappa=kappa,messages=messages)
	} else {
		res <- binomial.geo.Bayes.lr(formula=formula,units.m=units.m,coords=coords,
		           data=data,knots=knots,control.mcmc=control.mcmc,
		           control.prior=control.prior,kappa=kappa,messages=messages)		
	}
	res$call <- match.call()
	return(res)
}


##' @title Bayesian spatial prediction for the binomial logistic and binary probit models
##' @description This function performs Bayesian spatial prediction for the binomial logistic and binary probit models.
##' @param object an object of class "Bayes.PrevMap" obtained as result of a call to \code{\link{binomial.logistic.Bayes}} or \code{\link{binary.probit.Bayes}}.
##' @param grid.pred a matrix of prediction locations.
##' @param predictors a data frame of the values of the explanatory variables at each of the locations in \code{grid.pred}; each column correspond to a variable and each row to a location. \bold{Warning:} the names of the columns in the data frame must match those in the data used to fit the model. Default is \code{predictors=NULL} for models with only an intercept.
##' @param type a character indicating the type of spatial predictions: \code{type="marginal"} for marginal predictions or \code{type="joint"} for joint predictions. Default is \code{type="marginal"}. In the case of a low-rank approximation only joint predictions are available.
##' @param scale.predictions a character vector of maximum length 3, indicating the required scale on which spatial prediction is carried out: "logit", "prevalence", "odds" and "probit". Default is \code{scale.predictions="prevalence"}.
##' @param quantiles a vector of quantiles used to summarise the spatial predictions.
##' @param standard.errors logical; if \code{standard.errors=TRUE}, then standard errors for each \code{scale.predictions} are returned. Default is \code{standard.errors=FALSE}.
##' @param thresholds a vector of exceedance thresholds; default is \code{NULL}.
##' @param scale.thresholds a character value ("logit", "prevalence", "odds" or "probit") indicating the scale on which exceedance thresholds are provided.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return A "pred.PrevMap" object list with the following components: \code{logit}; \code{prevalence}; \code{odds}; \code{probit};\code{exceedance.prob}, corresponding to a matrix of the exceedance probabilities where each column corresponds to a specified value in \code{thresholds}; \code{samples}, corresponding to a matrix of the posterior samples at each prediction locations for the linear predictor; \code{grid.pred} prediction locations. 
##' Each of the three components \code{logit}, \code{prevalence},  \code{odds} and \code{probit} is also a list with the following components:
##' @return \code{predictions}: a vector of the predictive mean for the associated quantity (logit, odds or prevalence).
##' @return \code{standard.errors}: a vector of prediction standard errors (if \code{standard.errors=TRUE}).
##' @return \code{quantiles}: a matrix of quantiles of the resulting predictions with each column corresponding to a quantile specified through the argument \code{quantiles}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom geoR matern 
##' @export 
spatial.pred.binomial.Bayes <- function(object,grid.pred,predictors=NULL,
                               type="marginal",
                               scale.predictions="prevalence",
                               quantiles=c(0.025,0.975),
                               standard.errors=FALSE,thresholds=NULL,
                               scale.thresholds=NULL,messages=TRUE) {
    
    check.probit <- substr(object$call[1],8,13)=="probit"                         	
    if(length(scale.predictions) > 3) stop("too many values for scale.predictions")
    
    if(nrow(grid.pred) < 2) stop("prediction locations must be at least two.")
    if(length(predictors)>0 && class(predictors)!="data.frame") stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.") 
    if(length(predictors)>0 && any(is.na(predictors))) stop("missing values found in 'predictors'.")

    if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions should be marginal or joint")
	ck <- length(dim(object$knots)) > 0
	if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions must be marginal or joint.")
	if(ck & type=="marginal") warning("only joint predictions are available for the low-rank approximation")
    out <- list()
	for(i in 1:length(scale.predictions)) {
		if(any(c("logit","prevalence","odds","probit")==scale.predictions[i])==
		          FALSE) stop("invalid scale.predictions")
	}
	
	if(any(scale.predictions == "logit") & check.probit) {
		warning("a binary probit model was fitted, hence spatial prediction will be carried out on the probit scale and not on the logit scale")
	}
	if(length(thresholds)>0) {
	   if(any(c("logit","prevalence","odds","probit")
		  ==scale.thresholds)==FALSE) {
	      stop("scale thresholds must be logit, prevalence or odds scale.")
	   }
	}
    
    if(length(thresholds)==0 & length(scale.thresholds)>0 |
       length(thresholds)>0 & length(scale.thresholds)==0) stop("to estimate exceedance probabilities both thresholds and scale.thresholds must be provided.")
	p <- ncol(object$D)
	kappa <- object$kappa	
	n.pred <- nrow(grid.pred)
    D <- object$D
    
    
	if(p==1) {
	   predictors <- matrix(1,nrow=n.pred)	
	} else {
	   if(length(dim(predictors))==0) stop("covariates at prediction locations should be provided.")	
	   predictors <- as.matrix(model.matrix(
	                        delete.response(terms(formula(object$call))),data=predictors))
	   if(nrow(predictors)!=nrow(grid.pred)) stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")        
	   if(ncol(predictors)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")             
	}				
	
	n.samples <- nrow(object$estimate)	
	out <- list()
	eta.sim <- matrix(NA,nrow=n.samples,ncol=n.pred)
	fixed.nugget <- ncol(object$estimate) < p+3
	if(ck) {
	   	
	   U.k <- as.matrix(pdist(object$coords,object$knots))
	   U.k.pred <- as.matrix(pdist(grid.pred,object$knots))	
	   if(messages) cat("Bayesian prediction (this might be time consuming) \n")
	   for(j in 1:n.samples) {		
	   	   rho <- 2*sqrt(kappa)*object$estimate[j,"phi"]
		   mu.pred <- as.numeric(predictors%*%object$estimate[j,1:p])
		   K.pred <- matern.kernel(U.k.pred,rho,object$kappa)
		   eta.sim[j,] <- as.numeric(mu.pred+K.pred%*%object$S[j,])
		   if(messages) cat("Iteration ",j," out of ",n.samples,"\r")
	   }
	   if(messages) cat("\n")    
       if(any(scale.predictions=="logit") | 
          any(scale.predictions=="probit")) {
       	 if(messages) {
       	 	if(check.probit) {
       	 		cat("Spatial prediction: probit \n")
       	 	} else {
       	 		cat("Spatial prediction: logit \n")
       	 	}	
       	 }	
    	     if(check.probit) {
    	        out$probit$predictions <- apply(eta.sim,2,mean)
    	        if(length(quantiles)>0) out$probit$quantiles <- 
    	        t(apply(eta.sim,2,function(c) quantile(c,quantiles)))
    	        if(standard.errors) out$probit$standard.errors <- 
    	        apply(eta.sim,2,sd)	
    	     } else {
    	     	out$logit$predictions <- apply(eta.sim,2,mean)	
    	     	if(length(quantiles)>0) out$logit$quantiles <- 
    	        t(apply(eta.sim,2,function(c) quantile(c,quantiles)))
    	        if(standard.errors) out$logit$standard.errors <- 
    	        apply(eta.sim,2,sd)
    	     }
        }
    
        if(any(scale.predictions=="odds")) {
         if(messages) cat("Spatial prediction: odds \n")
    	     odds <- exp(eta.sim)
    	     out$odds$predictions <- apply(odds,2,mean)
    	     
    	     if(length(quantiles)>0) out$odds$quantiles <- 
    	     t(apply(odds,2,function(c) quantile(c,quantiles)))
    	     
    	     if(standard.errors) out$odds$standard.errors <-  
    	     apply(odds,2,sd)
        }              
          
       if(any(scale.predictions=="prevalence")) {
       	 if(messages) cat("Spatial prediction: prevalence \n")
    	     if(check.probit) {
    	     	prev <- exp(eta.sim)/(1+exp(eta.sim))
    	     } else {
    	     	prev <- pnorm(eta.sim)
    	     }
    	     out$prevalence$predictions <- apply(prev,2,mean)
    	     if(length(quantiles)>0) out$prevalence$quantiles <- 
    	     t(apply(prev,2,function(c) quantile(c,quantiles)))
    	     
    	     if(standard.errors) out$prevalence$standard.errors <- 
    	     apply(prev,2,sd)
        }
        
        if(length(scale.thresholds) > 0) {
        	 if(messages) cat("Spatial prediction: exceedance probabilities \n")
             if(scale.thresholds=="prevalence") {
                thresholds <- log(thresholds/(1-thresholds))	
             } else if(scale.thresholds=="odds") {
             	thresholds <- log(thresholds)
             } else if(scale.thresholds=="probit") {
             	thresholds <- qnorm(thresholds)
             }
             out$exceedance.prob <- sapply(thresholds, function(x) 
             apply(eta.sim,2,function(c) mean(c > x)))	
       } 
	} else {
	   U <- dist(object$coords)
       mu.cond <- matrix(NA,nrow=n.samples,ncol=n.pred)
       sd.cond <- matrix(NA,nrow=n.samples,ncol=n.pred)
       U.pred.coords <- as.matrix(pdist(grid.pred,object$coords))
	   if(type=="joint") {
	   	if(messages) cat("Type of predictions: joint \n")
	   	U.grid.pred <- dist(grid.pred)
	   } else {
	   	if(messages) cat("Type of predictions: marginal \n")
	   }

       for(j in 1:n.samples)  { 
          sigma2 <- object$estimate[j,"sigma^2"]
          phi <-object$estimate[j,"phi"]
          if(!fixed.nugget) {
             tau2 <- object$estimate[j,"tau^2"]
          } else {
       	     tau2 <- 0
          }
          mu <- D%*%object$estimate[j,1:p]    
   	      Sigma <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                   cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov 
	      Sigma.inv <- solve(Sigma)   
          C <- sigma2*matern(U.pred.coords,phi,kappa)	                
          A <- C%*%Sigma.inv
	      mu.pred <- as.numeric(predictors%*%object$estimate[j,1:p])         	             	
          if(length(object$ID.coords)>0) {
	         mu.cond[j,] <- mu.pred+A%*%object$S[j,]        
          } else { 	            
             mu.cond[j,] <- mu.pred+A%*%(object$S[j,]-mu)   	
          }
	              
          if(type=="marginal") {
      	     sd.cond[j,] <- sqrt(sigma2-diag(A%*%t(C)))
    	     eta.sim[j,] <-  rnorm(n.pred,mu.cond[j,],sd.cond[j,])
          } else if (type=="joint") {
    	     Sigma.pred <-  varcov.spatial(dists.lowertri=U.grid.pred,cov.model="matern",
	                                  cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
	         Sigma.cond <- Sigma.pred - A%*%t(C)
	         Sigma.cond.sroot <- t(chol(Sigma.cond))
	         sd.cond <- sqrt(diag(Sigma.cond))	   
	         eta.sim[j,] <- mu.cond[j,]+Sigma.cond.sroot%*%rnorm(n.pred)
         }  

         if(messages) cat("Iteration ",j," out of ",n.samples,"\r")   	             
      }
      cat("\n")    
      if(any(scale.predictions=="logit") | 
         any(scale.predictions=="probit")) {
      	if(messages) {
       	 	if(check.probit) {
       	 		cat("Spatial prediction: probit \n")
       	 	} else {
       	 		cat("Spatial prediction: logit \n")
       	 	}	
       	 }
    	if(check.probit) {
    	   out$probit$predictions <- apply(mu.cond,2,mean)
    	   
    	   if(length(quantiles)>0) out$probit$quantiles <- 
    	   t(apply(eta.sim,2,function(c) quantile(c,quantiles)))
    	   
    	   if(standard.errors) out$probit$standard.errors <-  
    	   apply(eta.sim,2,sd)	
    	} else {
           out$logit$predictions <- apply(mu.cond,2,mean)
    	   
    	   if(length(quantiles)>0) out$logit$quantiles <- 
    	   t(apply(eta.sim,2,function(c) quantile(c,quantiles)))
    	   
    	   if(standard.errors) out$logit$standard.errors <-  
    	   apply(eta.sim,2,sd)		
    	}  
      }
    
      if(any(scale.predictions=="odds")) {
       	 if(messages) cat("Spatial prediction: odds \n")
    	     odds <- exp(eta.sim)
    	     out$odds$predictions <- 
    	     apply(exp(mu.cond+0.5*(sd.cond^2)),2,mean)
    	     
    	 if(length(quantiles)>0) out$odds$quantiles <- 
    	     t(apply(odds,2,function(c) quantile(c,quantiles)))
    	     
    	 if(standard.errors) out$odds$standard.errors <-  
    	     apply(odds,2,sd)
      }
    
      if(any(scale.predictions=="prevalence")) {
       	 if(messages) cat("Spatial prediction: prevalence \n")
    	 if(check.probit) {
    	 	prev <- pnorm(eta.sim)
    	 } else {
    	    prev <- exp(eta.sim)/(1+exp(eta.sim))	
    	 }
    	 
    	 out$prevalence$predictions <- apply(prev,2,mean)
    	 
    	 if(length(quantiles)>0) out$prevalence$quantiles <- 
    	 t(apply(prev,2,function(c) quantile(c,quantiles)))
    	 
    	 if(standard.errors) out$prevalence$standard.errors <-  
    	 apply(prev,2,sd)
      }
        
      if(length(scale.thresholds) > 0) {
         if(messages) cat("Spatial prediction: exceedance probabilities \n")
         if(scale.thresholds=="prevalence") {
                thresholds <- log(thresholds/(1-thresholds))	
         } else if(scale.thresholds=="odds") {
             	thresholds <- log(thresholds)
         } else if(scale.thresholds=="probit") {
             	thresholds <- qnorm(thresholds)
         }
         out$exceedance.prob <- sapply(thresholds, function(x) 
         apply(eta.sim,2,function(c) mean(c > x)))	
       }       	
   }
   out$grid.pred <- grid.pred
   out$samples <- eta.sim
   class(out) <- "pred.PrevMap"	
   out
}

##' @title Auxliary function for controlling profile log-likelihood in the linear Gaussian model
##' @description Auxiliary function used by \code{\link{loglik.linear.model}}. This function defines whether the profile-loglikelihood should be computed or evaluation of the likelihood is required by keeping the other parameters fixed. 
##' @param phi a vector of the different values that should be used in the likelihood evalutation for the scale parameter \code{phi}, or \code{NULL} if a single value is provided either as first argument in \code{start.par} (for profile likelihood maximization) or as fixed value in \code{fixed.phi}; default is \code{NULL}.
##' @param rel.nugget a vector of the different values that should be used in the likelihood evalutation for the relative variance of the nugget effect \code{nu2}, or \code{NULL} if a single value is provided either in \code{start.par} (for profile likelihood maximization) or as fixed value in \code{fixed.nu2}; default is \code{NULL}.
##' @param fixed.beta a vector for the fixed values of the regression coefficients \code{beta}, or \code{NULL} if profile log-likelihood is to be performed; default is \code{NULL}.
##' @param fixed.sigma2 value for the fixed variance of the Gaussian process \code{sigma2}, or \code{NULL} if profile log-likelihood is to be performed; default is \code{NULL}.
##' @param fixed.phi value for the fixed scale parameter \code{phi} in the Matern function, or \code{NULL} if profile log-likelihood is to be performed; default is \code{NULL}.
##' @param fixed.rel.nugget value for the fixed relative variance of the nugget effect; \code{fixed.rel.nugget=NULL} if profile log-likelihood is to be performed; default is \code{NULL}.
##' @seealso \code{\link{loglik.linear.model}}
##' @return A list with components named as the arguments. 
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
control.profile <- function(phi=NULL,rel.nugget=NULL,
                                         fixed.beta=NULL,fixed.sigma2=NULL,
                                         fixed.phi=NULL,
                                         fixed.rel.nugget=NULL) {
      	                                     	
   if(length(phi)==0 & length(rel.nugget)==0) stop("either phi or rel.nugget must be provided.")                                      	
   if(length(phi) > 0) {
     if(any(phi < 0)) stop("phi must be positive.")
   }	
   
   if(length(rel.nugget) > 0) {
     if(any(rel.nugget < 0)) stop("rel.nugget must be positive.")
   }	
   
   if(length(fixed.beta)!=0) {
   	  if(length(fixed.rel.nugget)==0 & length(fixed.phi)==0 & length(rel.nugget)==0 & length(phi)==0) stop("missing fixed value for phi or rel.nugget.")
   }   
   
   if(length(fixed.sigma2)==1) {
     if(fixed.sigma2 < 0) stop("fixed.sigma2 must be positive")
     if(length(fixed.rel.nugget)==0 & length(fixed.phi)==0 & length(rel.nugget)==0 & length(phi)==0) stop("missing fixed value for phi or rel.nugget.")
   }
      
   fixed.par.phi <- FALSE
   fixed.par.rel.nugget <- FALSE
   
   if((length(fixed.beta)>0&length(fixed.sigma2)>0)) {
     if(length(fixed.phi) > 0 && fixed.phi < 0) stop("fixed.phi must be positive")
     if(length(fixed.phi)>0 & length(fixed.rel.nugget)==0) fixed.par.rel.nugget <- TRUE
     if(any(c(length(fixed.beta),length(fixed.sigma2))==0)) stop("fixed.beta and fixed.sigma2 must be provided for evaluation of the likelihood with fixed paramaters.")      
     if(length(fixed.rel.nugget) > 0 && fixed.rel.nugget < 0) stop("fixed.rel.nugget must be positive")
     if(length(fixed.phi)==0 & length(fixed.rel.nugget)>0) fixed.par.phi <- TRUE
     if(any(c(length(fixed.beta),length(fixed.sigma2))==0)) stop("fixed.beta and fixed.sigma2 must be provided for evaluation of the likelihood with fixed paramaters.")
   } 

   
   if(length(rel.nugget) > 0 & length(fixed.rel.nugget) > 0) stop("rel.nugget and fixed.rel.nugget cannot be both provided.")	
      if(length(phi) > 0 & length(fixed.phi) > 0) stop("phi and fixed.phi cannot be both provided.")	
   if(fixed.par.phi | fixed.par.rel.nugget) {
      cat("Control profile: parameters have been set for likelihood evaluation with fixed parameters. \n")
   } else if(fixed.par.phi==FALSE & fixed.par.rel.nugget==FALSE) {
      cat("Control profile: parameters have been set for evaluation of the profile log-likelihood. \n")
   }   
   out <- list(phi=phi,rel.nugget=rel.nugget,
                   fixed.phi=fixed.phi,fixed.rel.nugget=fixed.rel.nugget,
                   fixed.beta=fixed.beta,fixed.sigma2=fixed.sigma2,
                   fixed.par.phi=fixed.par.phi,
                   fixed.par.rel.nugget=fixed.par.rel.nugget)  
   return(out)                                    	
}

##' @title Profile log-likelihood or fixed parameters likelihood evaluation for the covariance parameters in the geostatistical linear model 
##' @description Computes profile log-likelihood, or evaluatesx likelihood keeping the other paramaters fixed, for the scale parameter \code{phi} of the Matern function and the relative variance of the nugget effect \code{nu2} in the linear Gaussian model.
##' @param object an object of class 'PrevMap', which is the fitted linear model obtained with the function \code{\link{linear.model.MLE}}.
##' @param control.profile control parameters obtained with \code{\link{control.profile}}.
##' @param plot.profile logical; if \code{TRUE} a plot of the computed profile likelihood is displayed.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return an object of class "profile.PrevMap" which is a list with the following values
##' @return \code{eval.points.phi}: vector of the values used for \code{phi} in the evaluation of the likelihood.
##' @return \code{eval.points.rel.nugget}: vector of the values used for \code{nu2} in the evaluation of the likelihood.
##' @return \code{profile.phi}: vector of the values of the likelihood function evaluated at \code{eval.points.phi}.
##' @return \code{profile.rel.nugget}: vector of the values of the likelihood function evaluated at \code{eval.points.rel.nugget}.
##' @return \code{profile.phi.rel.nugget}: matrix of the values of the likelihood function evaluated at \code{eval.points.phi} and \code{eval.points.rel.nugget}.
##' @return \code{fixed.par}: logical value; \code{TRUE} is the evaluation if the likelihood is carried out by fixing the other parameters, and \code{FALSE} if the computation of the profile-likelihood was performed instead.
##' @importFrom geoR varcov.spatial
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
loglik.linear.model <- function(object,control.profile,plot.profile=TRUE,
                                  messages=TRUE) {
    if(class(object)!="PrevMap") stop("object must be of class PrevMap.")
	y <- object$y
	n <- length(y)
	D <- object$D	
	kappa <- object$kappa
	if(length(control.profile$fixed.beta)>0 && length(control.profile$fixed.beta)!=ncol(D)) stop("invalid value for fixed.beta") 
	U <- dist(object$coords)
	if(length(object$fixed.nugget)>0) stop("the nugget effect must not be fixed when using geo.linear.MLE")
	if(control.profile$fixed.par.phi==FALSE & control.profile$fixed.par.rel.nugget==FALSE
	   & length(control.profile$phi) > 0) {
	   	start.par <- object$estimate["log(nu^2)"]
	   } else if(control.profile$fixed.par.phi==FALSE & control.profile$fixed.par.rel.nugget==FALSE
	   & length(control.profile$rel.nugget) > 0) {
	   	start.par <- object$estimate["log(phi)"]   
	}
	   
	if(control.profile$fixed.par.phi | control.profile$fixed.par.rel.nugget) {
	log.lik <- function(beta,sigma2,phi,kappa,nu2) {
	      V <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),kappa=kappa,
	              nugget=nu2)$varcov
	      V.inv <- solve(V)
	      diff.beta <- y-D%*%beta
	      q.f <- as.numeric(t(diff.beta)%*%V.inv%*%diff.beta)
	      ldet.V <- determinant(V)$modulus
	      as.numeric(-0.5*(n*log(sigma2)+ldet.V+q.f/sigma2))
    }		
    out <- list()
    if(control.profile$fixed.par.phi & control.profile$fixed.par.rel.nugget==FALSE) {
       	out$eval.points.phi <- control.profile$phi
       	out$profile.phi <- rep(NA,length(control.profile$phi))
      	for(i in 1:length(control.profile$phi)) {
    		   	flush.console()
    		   	cat("Evaluation for phi=",control.profile$phi[i],"\r",sep="")
    		    out$profile.phi[i] <- log.lik(control.profile$fixed.beta,
    	                                       control.profile$fixed.sigma2,control.profile$phi[i],
    	                                       kappa,control.profile$fixed.rel.nugget)    			                        
        }
        flush.console()
        if(plot.profile) {
            plot(out$eval.points.phi,out$profile.phi,type="l",
            main=expression("Log-likelihood for"~phi~"and remaining parameters fixed", ),
            xlab=expression(phi),ylab="log-likelihood")
        
        }
    } else if(control.profile$fixed.par.phi==FALSE & control.profile$fixed.par.rel.nugget) {
       	out$eval.points.rel.nugget <- control.profile$rel.nugget
       	out$profile.rel.nugget <- rep(NA,length(control.profile$rel.nugget))
      	for(i in 1:length(control.profile$rel.nugget)) {
    		   flush.console()  
    		   cat("Evaluation for rel.nugget=",control.profile$rel.nugget[i],"\r",sep="")
    		   out$profile.rel.nugget[i] <- log.lik(control.profile$fixed.beta,
    			                                  control.profile$fixed.sigma2,control.profile$fixed.phi,
    			                                  kappa,control.profile$rel.nugget[i])   			                        
        }
        flush.console() 
        if(plot.profile) {
              plot(out$eval.points.rel.nugget,out$profile.rel.nugget,type="l",
              main=expression("Log-likelihood for"~nu^2~"and remaining parameters fixed", ),
              xlab=expression(nu^2),ylab="log-likelihood")
        }
    } else if(control.profile$fixed.par.phi &
                 control.profile$fixed.par.rel.nugget) {
       	out$eval.points.phi <- control.profile$phi
       	out$eval.points.rel.nugget <- control.profile$rel.nugget
       	out$profile.phi.rel.nugget <- matrix(NA,length(control.profile$phi),length(control.profile$rel.nugget))
      	for(i in 1:length(control.profile$phi)) {
      		for(j in 1:length(control.profile$rel.nugget)) {
    		   	 flush.console()  
    		   	 cat("Evaluation for phi=",control.profile$phi[i],
    		             " and rel.nugget=",control.profile$rel.nugget[j],"\r",sep="")
    			 out$profile.phi.rel.nugget[i,j] <- log.lik(control.profile$fixed.beta,
    			                                  control.profile$fixed.sigma2,control.profile$phi[i],
    			                                  kappa,control.profile$rel.nugget[j])  	    			                        
    	  	}      	  	     		
        }
        flush.console()  
        if(plot.profile) {
               contour(out$eval.points.phi,out$eval.points.rel.nugget,
               out$profile.phi.rel.nugget,
               main=expression("Log-likelihood for"~phi~"and"~nu^2~
                            "and remaining parameters fixed"),
               xlab=expression(phi),ylab=expression(nu^2))	
        }
    }	
    } else  {			
    log.lik.profile <- function(phi,kappa,nu2) {
	      V <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),kappa=kappa,
	              nugget=nu2)$varcov
	      V.inv <- solve(V)
	      beta.hat <- as.numeric(solve(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%y)
	      diff.beta.hat <- y-D%*%beta.hat
	      sigma2.hat <- as.numeric(t(diff.beta.hat)%*%V.inv%*%diff.beta.hat/n)
	      ldet.V <- determinant(V)$modulus
	      as.numeric(-0.5*(n*log(sigma2.hat)+ldet.V))
    }		
    out <- list()
    
    if(length(control.profile$phi) > 0 & length(control.profile$rel.nugget)==0) {
       	out$eval.points.phi <- control.profile$phi
       	out$profile.phi <- rep(NA,length(control.profile$phi))
      	for(i in 1:length(control.profile$phi)) {
    		   	flush.console()  
    		   	cat("Evaluation for phi=",control.profile$phi[i],"\r",sep="")
    			out$profile.phi[i] <- -nlminb(start.par,function(x) -log.lik.profile(control.profile$phi[i],
    			                        kappa,exp(x)))$objective    			                        
        }
        flush.console()
        if(plot.profile) {
           plot(out$eval.points.phi,out$profile.phi,type="l",
                  main=expression("Profile log-likelihood for"~phi),
                  xlab=expression(log(phi)),ylab="log-likelihood")	
        }
    } else if(length(control.profile$phi) == 0 &
       length(control.profile$rel.nugget) > 0) {
       	out$eval.points.rel.nugget <- control.profile$rel.nugget
       	out$profile.rel.nugget <- rep(NA,length(control.profile$rel.nugget))
      	for(i in 1:length(control.profile$rel.nugget)) {
    		   	flush.console()  
    		   	cat("Evaluation for rel.nugget=",control.profile$rel.nugget[i],"\r",sep="")
    			out$profile.rel.nugget[i] <- -nlminb(start.par,
    	                       function(x) -log.lik.profile(exp(x),kappa,
                              control.profile$rel.nugget[i]))$objective   			                        
        }
        flush.console()
        if(plot.profile) {
               plot(out$eval.points.rel.nugget,out$profile.rel.nugget,type="l",
               main=expression("Profile log-likelihood for"~nu^2),
               xlab=expression(nu^2),ylab="log-likelihood")
        }
    } else if(length(control.profile$phi) > 0 &
       length(control.profile$rel.nugget) > 0) {
       	out$eval.points.phi <- control.profile$phi
       	out$eval.points.rel.nugget <- control.profile$rel.nugget
       	out$profile.phi.rel.nugget <- matrix(NA,length(control.profile$phi),length(control.profile$rel.nugget))
      	for(i in 1:length(control.profile$phi)) {
      		for(j in 1:length(control.profile$rel.nugget)) {
    		   	flush.console()  
    		   	cat("Evaluation for phi=",control.profile$phi[i],
    		             " and rel.nugget=",control.profile$rel.nugget[j],"\r",sep="")
    			out$profile.phi.rel.nugget[i,j] <- log.lik.profile(control.profile$phi[i],
    		                                            kappa,
    			                                        control.profile$rel.nugget[j])   	    			                        
    	  	}      	  	     		
        }
        if(plot.profile) {
               contour(out$eval.points.phi,out$eval.points.rel.nugget,
               out$profile.phi.rel.nugget,
               main=expression("Profile log-likelihood for"~phi~"and"~nu^2),
               xlab=expression(phi),ylab=expression(nu^2))	
        }
      }           
    }        
    out$fixed.par <- control.profile$fixed.par.phi | control.profile$fixed.par.rel.nugget
    class(out) <- "profile.PrevMap"            
    return(out)   
}

##' @title Plot of the profile log-likelihood for the covariance parameters of the Matern function
##' @description This function displays a plot of the profile log-likelihood that is computed by the function \code{\link{loglik.linear.model}}.
##' @param x object of class "profile.PrevMap" obtained as output from \code{\link{loglik.linear.model}}.
##' @param log.scale logical; if \code{log.scale=TRUE}, the profile likleihood is plotted on the log-scale of the parameter values.
##' @param plot.spline.profile logical; if \code{TRUE} an interpolating spline of the profile-likelihood of for a univariate parameter is plotted. Default is \code{FALSE}. 
##' @param ... further arugments passed to \code{\link{plot}} if the profile log-likelihood is for only one parameter, or to \code{\link{contour}} for the bi-variate profile-likelihood.
##' @return A plot is returned. No value is returned.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method plot profile.PrevMap
##' @export
plot.profile.PrevMap <- function(x,log.scale=FALSE,plot.spline.profile=FALSE,...) {
	if(class(x) != "profile.PrevMap") stop("x must be of class profile.PrevMap")
    if(class(x$profile)=="matrix" & plot.spline.profile==TRUE) warning("spline interpolation is available only for univariate profile-likelihood")    
	if(class(x$profile)=="numeric") {
		plot.list <- list(...)
		if(length(plot.list$type)==0) plot.list$type <- "l"
		if(length(plot.list$xlab)==0) plot.list$xlab <- ""
		if(length(plot.list$ylab)==0) plot.list$ylab <- ""
        plot.list$x <- x[[1]]
        plot.list$y <- x[[2]]
		if(log.scale) {		
		   plot.list$x <- log(plot.list$x)	
		   do.call(plot,plot.list)
		   if(plot.spline.profile) {
		      f <- splinefun(log(x[[1]]),x[[2]])
			  lines(log(x[[1]]),f,col=2)	
		   }
		} else {
		   do.call(plot,plot.list)
		   if(plot.spline.profile) {
		      f <- splinefun(x[[1]],x[[2]])
			  lines(x[[1]],f,col=2)	
		   }
		}		    					
   	} else if(class(x$profile)=="matrix") {
   	  	if(log.scale) {
   	  	    contour(log(x[[1]]),log(x[[2]]),x[[3]],...)
  	  	} else {
            contour(x[[1]],x[[2]],x[[3]],...)   	  		
   	  	}
	 }				
}

##' @title Profile likelihood confidence intervals
##' @description Computes confidence intervals based on the interpolated profile likelihood computed for a single covariance parameter.
##' @param object object of class "profile.PrevMap" obtained from \code{\link{loglik.linear.model}}.
##' @param coverage a value between 0 and 1 indicating the coverage of the confidence interval based on the interpolated profile likelihood. Default is \code{coverage=0.95}.
##' @param plot.spline.profile logical; if \code{TRUE} an interpolating spline of the profile-likelihood of for a univariate parameter is plotted. Default is \code{FALSE}. 
##' @return A list with elements \code{lower} and \code{upper} for the upper and lower limits of the confidence interval, respectively.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export

loglik.ci <- function(object,coverage=0.95,plot.spline.profile=TRUE) {
	if(class(object) != "profile.PrevMap") stop("object must be of class profile.PrevMap")
	if(coverage <= 0 | coverage >=1) stop("'coverage' must be between 0 and 1.")
	if(object$fixed.par | class(object$profile)=="matrix") {
		stop("profile likelihood must be computed and only for a single paramater.")
	} else {
		f <- splinefun(object[[1]],object[[2]])
        max.log.lik <- -optimize(function(x)-f(x),
                                lower=min(object[[1]]),upper=max(object[[1]]))$objective
        par.set <- seq(min(object[[1]]),max(object[[1]]),length=20)
        k <- qchisq(coverage,df=1)
        par.val <- 2*(max.log.lik-f(par.set))-k
        done <- FALSE      
        if(par.val[1] < 0 & par.val[length(par.val)] > 0) {
        	   stop("smaller value of the parameter should be included when computing the profile likelihood")
        } else if(par.val[1] > 0 & par.val[length(par.val)] < 0) {
        	   stop("larger value of the parameter should be included when computing the profile likelihood")
        } else if(par.val[1] < 0 & par.val[length(par.val)] < 0) {
           stop("both smaller and larger value of the parameter should be included when computing the profile likelihood")        	
        }
        
        i <- 1
        lower.int <- c(NA,NA)
        upper.int <- c(NA,NA)
        out <- list()
        while(!done) {
        	    i <- i+1           	
           	if(sign(par.val[i-1]) != sign(par.val[i]) & sign(par.val[i-1])==1) {
           		lower.int[1] <- par.set[i-1]
           		lower.int[2] <- par.set[i+1]
           	}           	
           	if(sign(par.val[i-1]) != sign(par.val[i]) & sign(par.val[i-1])==-1) {
           		upper.int[1] <- par.set[i-1]
           		upper.int[2] <- par.set[i+1]
           		done <- TRUE
           	}       
        }
        out$lower <- uniroot(function(x) 2*(max.log.lik-f(x))-k,lower=lower.int[1],upper=lower.int[2])$root                        
        out$upper <- uniroot(function(x) 2*(max.log.lik-f(x))-k,lower=upper.int[1],upper=upper.int[2])$root 
        if(plot.spline.profile) {    
           a <- seq(min(par.set),max(par.set),length=1000)  
           b <- sapply(a, function(x) -2*(max.log.lik-f(x)))  
           plot(a,b,type="l",ylab="2*(profile log-likelihood - max log-likelihood)",
           xlab="parameter values",
           main="")
           abline(h=-k,lty="dashed")
           abline(v=out$lower,lty="dashed")
           abline(v=out$upper,lty="dashed")
        }
        cat("Likelihood-based ",100*coverage,"% confidence interval: (",out$lower,", ",out$upper,") \n",sep="")
	}
	out
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom geoR varcov.spatial
geo.linear.Bayes <- function(formula,coords,data,
                                              control.prior,
                                              control.mcmc,
                                              kappa,messages) {
    if(any(is.na(data))) stop("missing values are not accepted")                                                                                             
   	coords <- as.matrix(model.frame(coords,data))
	if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")
   kappa <- as.numeric(kappa)
   if(length(control.prior$log.prior.nugget)!=length(control.mcmc$start.nugget)) {
      stop("missing prior or missing starting value for the nugget effect")	
   }
   
   out <- list()

   mf <- model.frame(formula,data=data)
   y <- as.numeric(model.response(mf))
   n <- length(y)
   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))   
   p <- ncol(D)
   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   U <- dist(coords)
   V.inv <- control.prior$V.inv
   m <- control.prior$m   	  
   nu <- 2*kappa
   log.prior.theta1_2 <- function(theta1,theta2) {
         sigma2 <- exp(2*theta1)
         phi <- (exp(2*theta1-theta2))^(1/nu)
         log(phi)+log(sigma2)+control.prior$log.prior.sigma2(sigma2)+
         control.prior$log.prior.phi(phi)	
   }
      
   log.prior.theta3 <- function(theta3) {
      	 tau2 <- exp(theta3)
      	 theta3+control.prior$log.prior.nugget(tau2)
   }
      
   fixed.nugget <- length(control.mcmc$start.nugget)==0
       
   log.posterior <- function(theta1,theta2,beta,param.post,theta3=NULL) {
	     sigma2 <- exp(2*theta1)
	     mu <- D%*%beta
	     diff.beta <- beta-m
	     if(length(theta3)==1) {
	     	lp.theta3 <- log.prior.theta3(theta3)
	     } else {
	     	lp.theta3 <- 0
	     }	
	     diff.y <- as.numeric(y-mu)
	     out <- log.prior.theta1_2(theta1,theta2)+lp.theta3+
	     -0.5*(p*log(sigma2)+t(diff.beta)%*%V.inv%*%diff.beta/sigma2)+
	     -0.5*(n*log(sigma2)+param.post$ldetR+
	     t(diff.y)%*%param.post$R.inv%*%(diff.y)/sigma2)    
	     as.numeric(out)     
   }
      
   acc.theta1 <- 0
   acc.theta2 <- 0
   if(fixed.nugget==FALSE) acc.theta3 <- 0
   n.sim <- control.mcmc$n.sim 
   burnin <- control.mcmc$burnin
   thin <- control.mcmc$thin
   h.theta1 <- control.mcmc$h.theta1
   h.theta2 <- control.mcmc$h.theta2
   h.theta3 <- control.mcmc$h.theta3
   
   c1.h.theta1 <- control.mcmc$c1.h.theta1
   c2.h.theta1 <- control.mcmc$c2.h.theta1
   c1.h.theta2 <- control.mcmc$c1.h.theta2
   c2.h.theta2 <- control.mcmc$c2.h.theta2
   c1.h.theta3 <- control.mcmc$c1.h.theta3
   c2.h.theta3 <- control.mcmc$c2.h.theta3
   if(!fixed.nugget) {
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+3)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi","tau^2")               
         tau2.curr <- control.mcmc$start.nugget
         theta3.curr <- log(tau2.curr)
    } else {  
      	 out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
         colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")    
         theta3.curr <- NULL
	  	 tau2.curr <- 0
    }    
      
    # Initialise all the parameters
    sigma2.curr <- control.mcmc$start.sigma2
    phi.curr <- control.mcmc$start.phi	  
    theta1.curr <- 0.5*log(sigma2.curr)
    theta2.curr <- log(sigma2.curr/(phi.curr^nu))
    beta.curr <- control.mcmc$start.beta
	
    # Compute the log-posterior density
	R.curr <- varcov.spatial(dists.lowertri=U,cov.model="matern",
                                          cov.pars=c(1,phi.curr),
                                          kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov

	param.post.curr <- list()
	param.post.curr$R.inv <- solve(R.curr)
	param.post.curr$ldetR <- determinant(R.curr)$modulus
	lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                    param.post.curr,theta3.curr)
                    
    h1 <- rep(NA,n.sim)
    h2 <- rep(NA,n.sim)
    if(!fixed.nugget) h3 <- rep(NA,n.sim)                    
    for(i in 1:n.sim) {
      	
        # Update theta1 
    	theta1.prop <- theta1.curr+h.theta1*rnorm(1)
    	sigma2.prop <- exp(2*theta1.prop)
	    phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
    	R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.prop)$varcov 
         
        param.post.prop <- list()	     
	    param.post.prop$R.inv <- solve(R.prop)
	    param.post.prop$ldetR <- determinant(R.prop)$modulus
	    lp.prop <- log.posterior(theta1.prop,theta2.curr,beta.curr,
	                                          param.post.prop,theta3.curr) 	                                         	     
	    if(log(runif(1)) < lp.prop-lp.curr) {
	       theta1.curr <- theta1.prop
	       param.post.curr <- param.post.prop
	       lp.curr <- lp.prop
	       sigma2.curr <- sigma2.prop
	       phi.curr <- phi.prop
	       acc.theta1 <- acc.theta1+1
	    } 
	    rm(theta1.prop,sigma2.prop,phi.prop)   	                                          
        h1[i] <- h.theta1 <- max(0,h.theta1 + 
	                       (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                 
    	     
    	# Update theta2 
    	theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	    phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
    	R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                         cov.pars=c(1,phi.prop),
	                                         kappa=kappa,nugget=tau2.curr/sigma2.curr)$varcov 
         
        param.post.prop <- list()	     
	    param.post.prop$R.inv <- solve(R.prop)
	    param.post.prop$ldetR <- determinant(R.prop)$modulus
	    lp.prop <- log.posterior(theta1.curr,theta2.prop,beta.curr,
	                                          param.post.prop,theta3.curr) 	                                         
	     
	    if(log(runif(1)) < lp.prop-lp.curr) {
	        theta2.curr <- theta2.prop
	        param.post.curr <- param.post.prop
	        lp.curr <- lp.prop
	        phi.curr <- phi.prop
	        acc.theta2 <- acc.theta2+1
	    }    	                      
	    rm(theta2.prop,phi.prop)                     
        h2[i] <- h.theta2 <- max(0,h.theta2 + 
	                       (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))  
	                                            
        # Update theta3
        if(!fixed.nugget) { 
            theta3.prop <- theta3.curr+h.theta3*rnorm(1)
	        tau2.prop <- exp(theta3.prop)
    	    R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                                        cov.pars=c(1,phi.curr),
	                                        kappa=kappa,nugget=tau2.prop/sigma2.curr)$varcov 
         
            param.post.prop <- list()	     
	        param.post.prop$R.inv <- solve(R.prop)
	        param.post.prop$ldetR <- determinant(R.prop)$modulus
	        lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                          param.post.prop,theta3.prop) 	                                         
	     
	        if(log(runif(1)) < lp.prop-lp.curr) {
	           theta3.curr <- theta3.prop
	           param.post.curr <- param.post.prop
	           lp.curr <- lp.prop
	           tau2.curr <- tau2.prop
	           acc.theta3 <- acc.theta3+1
	        }    	                      
	        rm(theta3.prop,tau2.prop)                     
            h3[i] <- h.theta3 <- max(0,h.theta3 + 
	                       (c1.h.theta3*i^(-c2.h.theta3))*(acc.theta3/i-0.45)) 
	     }               
	                          
         # Update beta
         A <- t(D)%*%param.post.curr$R.inv
         cov.beta <- solve(V.inv+A%*%D)
         mean.beta <- as.numeric(cov.beta%*%(V.inv%*%m+A%*%y))
         beta.curr <- as.numeric(mean.beta+sqrt(sigma2.curr)*
                            t(chol(cov.beta))%*%rnorm(p))
         lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                                 param.post.curr,theta3.curr)                         
      
         if(i > burnin & (i-burnin)%%thin==0) {
		    if(fixed.nugget) {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr)
		    	} else {
		    	   out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		    	                                                        phi.curr,tau2.curr)
		    }	
   	    }
        if(messages) cat("Iteration",i,"out of",n.sim,"\r")
        flush.console()                            
   }  	   
   class(out) <- "Bayes.PrevMap"	
   out$y <- y
   out$D <- D
   out$coords <- coords
   out$kappa <- kappa
   out$h1 <- h1
   out$h2 <- h2
   if(!fixed.nugget) out$h3 <- h3   
   return(out)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
geo.linear.Bayes.lr <- function(formula,coords,knots,data,
                                              control.prior,
                                              control.mcmc,
                                              kappa,messages) {
      knots <- as.matrix(knots)   
      if(any(is.na(data))) stop("missing values are not accepted")                                                                                                                             
      coords <- as.matrix(model.frame(coords,data))
	  if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")
      kappa <- as.numeric(kappa)
      if(length(control.prior$log.prior.nugget)!=
         length(control.mcmc$start.nugget)) {
         stop("missing prior or missing starting value for the nugget effect")	
      }
   
      out <- list()
      if(length(control.prior$log.prior.nugget)==0) stop("nugget effect must be included in the model.")
      mf <- model.frame(formula,data=data)
      y <- as.numeric(model.response(mf))
      n <- length(y)
      N <- nrow(knots)
      D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))   
      p <- ncol(D)
      if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
      U.k <- as.matrix(pdist(coords,knots))    

      V.inv <- control.prior$V.inv
      m <- control.prior$m   	  
      nu <- 2*kappa
      
      log.prior.theta1_2 <- function(theta1,theta2) {
           sigma2 <- exp(2*theta1)
           phi <- (exp(2*theta1-theta2))^(1/nu)
           log(phi)+log(sigma2)+control.prior$log.prior.sigma2(sigma2)+
           control.prior$log.prior.phi(phi)	
       }
      
       log.prior.theta3 <- function(theta3) {
      	   tau2 <- exp(theta3)
      	   theta3+control.prior$log.prior.nugget(tau2)
       }
             
       log.posterior <- function(theta1,theta2,beta,S,W,theta3) {
	       sigma2 <- exp(2*theta1)
	       tau2 <- exp(theta3)
	       mu <- D%*%beta
	       mean.y <- mu+W
	       diff.beta <- beta-m	    
	       lp.theta3 <- log.prior.theta3(theta3)
	       diff.y <- as.numeric(y-mean.y)	       
	       out <- log.prior.theta1_2(theta1,theta2)+lp.theta3+
	                 -0.5*(t(diff.beta)%*%V.inv%*%diff.beta)+
	                 -0.5*(N*log(sigma2)+sum(S^2)/sigma2)
	                 -0.5*(n*log(tau2)+sum(diff.y^2)/tau2)    
	       return(as.numeric(out))
       }
      
       acc.theta1 <- 0
       acc.theta2 <- 0
       acc.theta3 <- 0
       n.sim <- control.mcmc$n.sim 
       burnin <- control.mcmc$burnin
       thin <- control.mcmc$thin
       h.theta1 <- control.mcmc$h.theta1
       h.theta2 <- control.mcmc$h.theta2
       h.theta3 <- control.mcmc$h.theta3
       
       c1.h.theta1 <- control.mcmc$c1.h.theta1
       c2.h.theta1 <- control.mcmc$c2.h.theta1
       c1.h.theta2 <- control.mcmc$c1.h.theta2
       c2.h.theta2 <- control.mcmc$c2.h.theta2
       c1.h.theta3 <- control.mcmc$c1.h.theta3
       c2.h.theta3 <- control.mcmc$c2.h.theta3
       out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+3)	
       colnames(out$estimate) <- c(colnames(D),"sigma^2","phi","tau^2")
       out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=N)
       out$const.sigma2 <- rep(NA,(n.sim-burnin)/thin)             
      
       # Initialise all the parameters
       beta.curr <- control.mcmc$start.beta
       sigma2.curr <- control.mcmc$start.sigma2
	   phi.curr <- control.mcmc$start.phi
	   rho.curr <- phi.curr*2*sqrt(kappa)
       tau2.curr <- control.mcmc$start.nugget
       theta3.curr <- log(tau2.curr)
	   nu2.curr <- tau2.curr/sigma2.curr
	   theta1.curr <- 0.5*log(sigma2.curr)
	   theta2.curr <- log(sigma2.curr/(phi.curr^nu))
	   S.curr <- rep(0,N)
	   W.curr <- rep(0,n)
	
       # Compute the log-posterior density
  	   K.curr <- matern.kernel(U.k,rho.curr,kappa)
	   lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                       S.curr,W.curr,theta3.curr)
       h1 <- rep(NA,n.sim)
       h2 <- rep(NA,n.sim)
       h3 <- rep(NA,n.sim)    
       for(i in 1:n.sim) {             	 	
          
           # Update theta1 
    	   theta1.prop <- theta1.curr+h.theta1*rnorm(1)
    	   sigma2.prop <- exp(2*theta1.prop)
	       phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
           rho.prop <- phi.prop*2*sqrt(kappa)
           K.prop <- matern.kernel(U.k,rho.prop,kappa)
           W.prop <- as.numeric(K.prop%*%S.curr)
	       lp.prop <- log.posterior(theta1.prop,theta2.curr,beta.curr,
	                                            S.curr,W.prop,theta3.curr) 	                                         	     
	       if(log(runif(1)) < lp.prop-lp.curr) {
	          theta1.curr <- theta1.prop
	          W.curr <- W.prop
	          lp.curr <- lp.prop
	          sigma2.curr <- sigma2.prop
	          phi.curr <- phi.prop
	          rho.curr <- rho.prop
	          acc.theta1 <- acc.theta1+1
	          K.curr <- K.prop
	       } 
	       rm(theta1.prop,sigma2.prop,phi.prop,rho.prop,K.prop)   	                                          
           h1[i] <- h.theta1 <- max(0,h.theta1 + 
	                       (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                    
    	     
       	   # Update theta2 
    	   theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	       phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
           rho.prop <- phi.prop*2*sqrt(kappa)
           param.post.prop <- list()	     
           K.prop <- matern.kernel(U.k,rho.prop,kappa)
           W.prop <- as.numeric(K.prop%*%S.curr)           
	       lp.prop <- log.posterior(theta1.curr,theta2.prop,beta.curr,
	                                            S.curr,W.prop,theta3.curr) 	                                         
	     
	       if(log(runif(1)) < lp.prop-lp.curr) {
	           theta2.curr <- theta2.prop
	           W.curr <- W.prop
	           lp.curr <- lp.prop
	           phi.curr <- phi.prop
	           rho.curr <- rho.prop
	           acc.theta2 <- acc.theta2+1
	           K.curr <- K.prop
	       }    	                      
	       rm(theta2.prop,phi.prop,rho.prop)                     
           h2[i] <- h.theta2 <- max(0,h.theta2 + 
	                       (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))  
	                                            
           # Update theta3
           theta3.prop <- theta3.curr+h.theta3*rnorm(1)
	       tau2.prop <- exp(theta3.prop)
	       lp.prop <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                            S.curr,W.curr,theta3.prop) 	                                         
	     
	       if(log(runif(1)) < lp.prop-lp.curr) {
	           theta3.curr <- theta3.prop
	           lp.curr <- lp.prop
	           tau2.curr <- tau2.prop
	           acc.theta3 <- acc.theta3+1
	       }    	                      
	       rm(theta3.prop,tau2.prop)                     
           h3[i] <- h.theta3 <- max(0,h.theta3 + 
	                       (c1.h.theta3*i^(-c2.h.theta3))*(acc.theta3/i-0.45)) 
	                    
	                          
           # Update beta
           cov.beta <- solve(V.inv+t(D)%*%D/tau2.curr)
           mean.beta <- as.numeric(cov.beta%*%
                                 (V.inv%*%m+t(D)%*%(y-W.curr)/tau2.curr))
           beta.curr <- as.numeric(mean.beta+
                            t(chol(cov.beta))%*%rnorm(p))
           lp.curr <- log.posterior(theta1.curr,theta2.curr,beta.curr,
	                                           S.curr,W.curr,theta3.curr)                         
      
           # Update S
           cov.S <- t(K.curr)%*%K.curr/tau2.curr
           diag(cov.S) <- diag(cov.S)+1/sigma2.curr
           cov.S <- solve(cov.S)
           mean.S <- as.numeric(cov.S%*%t(K.curr)%*%
                             (y-D%*%beta.curr)/tau2.curr)
           S.curr <- as.numeric(mean.S+t(chol(cov.S))%*%rnorm(N))
           W.curr <- as.numeric(K.curr%*%S.curr)
           
           if(i > burnin & (i-burnin)%%thin==0) {
		    out$estimate[(i-burnin)/thin,]  <- c(beta.curr,sigma2.curr,phi.curr,tau2.curr)	
		    out$S[(i-burnin)/thin,] <- S.curr
		    out$const.sigma2[(i-burnin)/thin] <- mean(apply(
		                                                   K.curr,1,function(r) sqrt(sum(r^2))))	
   	       }
           if(messages) cat("Iteration",i,"out of",n.sim,"\r")
           flush.console()                            
   }   
   class(out) <- "Bayes.PrevMap"	
   out$y <- y
   out$D <- D
   out$coords <- coords
   out$kappa <- kappa  
   out$knots <- knots 
   out$h1 <- h1
   out$h2 <- h2
   out$h3 <- h3           
   return(out)
}

##' @title Bayesian spatial predictions for the geostatistical Linear Gaussian model 
##' @description This function performs Bayesian prediction for a geostatistical linear Gaussian model.
##' @param object an object of class "Bayes.PrevMap" obtained as result of a call to \code{\link{linear.model.Bayes}}.
##' @param grid.pred a matrix of prediction locations.
##' @param predictors a data frame of the values of the explanatory variables at each of the locations in \code{grid.pred}; each column correspond to a variable and each row to a location. \bold{Warning:} the names of the columns in the data frame must match those in the data used to fit the model. Default is \code{predictors=NULL} for models with only an intercept.
##' @param type a character indicating the type of spatial predictions: \code{type="marginal"} for marginal predictions or \code{type="joint"} for joint predictions. Default is \code{type="marginal"}. In the case of a low-rank approximation only joint predictions are available.
##' @param scale.predictions a character vector of maximum length 3, indicating the required scale on which spatial prediction is carried out: "logit", "prevalence" and "odds". Default is \code{scale.predictions=c("logit","prevalence","odds")}.
##' @param quantiles a vector of quantiles used to summarise the spatial predictions.
##' @param standard.errors logical; if \code{standard.errors=TRUE}, then standard errors for each \code{scale.predictions} are returned. Default is \code{standard.errors=FALSE}.
##' @param thresholds a vector of exceedance thresholds; default is \code{thresholds=NULL}.
##' @param scale.thresholds a character value indicating the scale on which exceedance thresholds are provided: \code{"logit"}, \code{"prevalence"} or \code{"odds"}. Default is \code{scale.thresholds=NULL}.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return A "pred.PrevMap" object list with the following components: \code{logit}; \code{prevalence}; \code{odds}; \code{exceedance.prob}, corresponding to a matrix of the exceedance probabilities where each column corresponds to a specified value in \code{thresholds}; \code{grid.pred} prediction locations. 
##' Each of the three components \code{logit}, \code{prevalence} and  \code{odds} is also a list with the following components:
##' @return \code{predictions}: a vector of the predictive mean for the associated quantity (logit, odds or prevalence).
##' @return \code{standard.errors}: a vector of prediction standard errors (if \code{standard.errors=TRUE}).
##' @return \code{quantiles}: a matrix of quantiles of the resulting predictions with each column corresponding to a quantile specified through the argument \code{quantiles}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}  
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom geoR matern 
##' @export 
spatial.pred.linear.Bayes <- function(object,grid.pred,predictors=NULL,
                                                    type="marginal",
                                                    scale.predictions=c("logit",
                                                    "prevalence","odds"),
                                                    quantiles=c(0.025,0.975),
                                                    standard.errors=FALSE,
                                                    thresholds=NULL,scale.thresholds=NULL,
                                                    messages=TRUE) {
    if(nrow(grid.pred) < 2) stop("prediction locations must be at least two.")
    if(length(predictors)>0 && class(predictors)!="data.frame") stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.")
    if(length(predictors)>0 && any(is.na(predictors))) stop("missing values found in 'predictors'.")

	object$p <- ncol(object$D)
	object$fixed.nugget <- ncol(object$estimate) < object$p+3
	kappa <- object$kappa	
	n.pred <- nrow(grid.pred)
	coords <- object$coords
	ck <- length(dim(object$knots)) > 0
	if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions must be marginal or joint.")
	if(ck & type=="marginal") warning("only joint predictions are available for the low-rank approximation")
    out <- list()
	for(i in 1:length(scale.predictions)) {
		if(any(c("logit","prevalence","odds")==scale.predictions[i])==
		          FALSE) stop("invalid scale.predictions")
	}
	
	if(length(thresholds)>0) {
		if(any(c("logit","prevalence","odds")==scale.thresholds)==FALSE) {
			stop("scale thresholds must be logit, prevalence or odds scale.")
		}
	}
    
    if(length(thresholds)==0 & length(scale.thresholds)>0 |
       length(thresholds)>0 & length(scale.thresholds)==0) stop("to estimate exceedance probabilities both thresholds and scale.thresholds must be provided.")
	if(object$p==1) {
	   predictors <- matrix(1,nrow=n.pred)	
	} else {
	   if(length(dim(predictors))==0) stop("covariates at prediction locations should be provided.")	
	   predictors <- as.matrix(model.matrix(delete.response(terms(formula(object$call))),data=predictors))
	   if(nrow(predictors)!=nrow(grid.pred)) stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")
	   if(ncol(predictors)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
	}
        
   	n.samples <- nrow(object$estimate)
	n.pred <- nrow(grid.pred)
    
	if(ck) {	
	  knots <- object$knots	
      U.k.pred.coords <- as.matrix(pdist(grid.pred,knots))
	} else {	
	  U <- dist(coords)
	  U.pred.coords <- as.matrix(pdist(grid.pred,coords))
    }
    
	predictions <- matrix(NA,n.samples,ncol=n.pred)	
	mu.cond <- matrix(NA,n.samples,ncol=n.pred)
    sd.cond <- matrix(NA,n.samples,ncol=n.pred)
 
	for(j in 1:n.samples) {	  	
	   beta <- object$estimate[j,1:object$p]
	   sigma2 <- object$estimate[j,"sigma^2"]
	   phi <- object$estimate[j,"phi"]
	   if(object$fixed.nugget){		
	        tau2 <- 0
   	    } else {
	        tau2 <- object$estimate[j,"tau^2"]
	   }
	   mu.pred <- as.numeric(predictors%*%beta)  
	   if(ck) {
	   	  mu <- as.numeric(object$D%*%beta)
	   	  rho <- 2*sqrt(kappa)*phi
	   	  K.pred <- matern.kernel(U.k.pred.coords,rho,kappa)
	   	  predictions[j,] <- mu.pred+K.pred%*%object$S[j,]
	   	  flush.console()
          if(messages) cat("Iteration ",j," out of ",n.samples," \r",sep="")
	   } else {
          Sigma <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                   cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov 
	      Sigma.inv <- solve(Sigma)  	
	      C <- sigma2*matern(U.pred.coords,phi,kappa)   		   
          A <- C%*%Sigma.inv
	      mu <- object$D%*%beta  	    		 	   	
       
          mu.cond[j,] <- as.numeric(mu.pred+A%*%(object$y-mu))    
          if(type=="marginal") {
               sd.cond[j,] <- sqrt(sigma2-diag(A%*%t(C)))
               predictions[j,] <- rnorm(n.pred,mu.cond[j,],sd.cond[j,]) 
	      } else if(type=="joint" & any(scale.predictions!="logit")) {
	   	       Sigma.pred <-  varcov.spatial(coords=grid.pred,cov.model="matern",
	           cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
	           Sigma.cond <- Sigma.pred - A%*%t(C)
	           Sigma.cond.sroot <- t(chol(Sigma.cond))
	           sd.cond[j,] <- sqrt(diag(Sigma.cond))
	           predictions[j,] <- mu.cond[j,]+Sigma.cond.sroot%*%rnorm(n.pred)
          }
          flush.console()
          if(messages) cat("Iteration ",j," out of ",n.samples," \r",sep="")
         }
       } 
       if(messages) cat("\n")               
       if(any(scale.predictions=="logit")) {
       	     if(messages) cat("Spatial prediction: logit \n")
       	     if(ck) {
       	     	out$logit$predictions <- apply(predictions,2,mean)
       	     } else {
       	        out$logit$predictions <- apply(mu.cond,2,mean)	
       	     }    	     
    	     if(length(quantiles)>0) out$logit$quantiles <- 
    	                                  t(apply(predictions,2,function(c) quantile(c,quantiles)))
    	     if(standard.errors) out$logit$standard.errors <-  apply(predictions,2,sd)
       }
    
       if(any(scale.predictions=="odds")) {
       	     if(messages) cat("Spatial prediction: odds \n")
    	     odds <- exp(predictions)
    	     if(ck) {
    	        out$odds$predictions <- apply(odds,2,mean)
    	     } else {
                out$odds$predictions <- apply(exp(mu.cond+0.5*(sd.cond^2)),2,mean)
    	     }
    	     if(length(quantiles)>0) out$odds$quantiles <- 
    	                                    t(apply(odds,2,function(c) quantile(c,quantiles)))
    	     if(standard.errors) out$odds$standard.errors <-  apply(odds,2,sd)
       }
          
       if(any(scale.predictions=="prevalence")) {
       	     if(messages) cat("Spatial prediction: prevalence \n")
    	     prev <- exp(predictions)/(1+exp(predictions))
    	     out$prevalence$predictions <- apply(prev,2,mean)
    	     if(length(quantiles)>0) out$prevalence$quantiles <- 
    	                                    t(apply(prev,2,function(c) quantile(c,quantiles)))
    	     if(standard.errors) out$prevalence$standard.errors <-  apply(prev,2,sd)
        }
        
        if(length(scale.thresholds) > 0) {
        	 if(messages) cat("Spatial prediction: exceedance probabilities \n")
             if(scale.thresholds=="prevalence") {
                thresholds <- log(thresholds/(1-thresholds))	
             } else if(scale.thresholds=="odds") {
             	thresholds <- log(thresholds)
             }
             out$exceedance.prob <- sapply(thresholds, function(x) 
                                                             apply(predictions,2,function(c) mean(c > x)))	
       }   
       out$grid.pred <- grid.pred
       out$samples <- predictions             	
       class(out) <- "pred.PrevMap"
       out
}


##' @title Extract model coefficients 
##' @description \code{coef} extracts parameters estimates from models fitted with the functions \code{\link{linear.model.MLE}} and \code{\link{binomial.logistic.MCML}}.
##' @param object an object of class "PrevMap".
##' @param ... other arguments.
##' @return coefficients extracted from the model object \code{object}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
coef.PrevMap <- function(object,...) {
   if(class(object)!="PrevMap") stop("object must be of class PrevMap.")	
   object$p <- ncol(object$D)
   out <- object$estimate
   if(length(dim(object$knots)) > 0 & length(object$units.m) > 0) {
      object$fixed.rel.nugget <- 0	
   }	
   out[-(1:object$p)] <- exp(out[-(1:object$p)])
   names(out)[object$p+1] <- "sigma^2"
   names(out)[object$p+2] <- "phi"	
	
   linear.ID.coords <- length(object$ID.coords) > 0 & 
                       substr(object$call[1],1,6)=="linear"

   if(linear.ID.coords & length(object$fixed.rel.nugget)==1) {
      out[object$p+3] <- out[object$p+1]*out[object$p+3]
      names(out)[object$p+3] <- "omega^2"
   }

   if(length(object$fixed.rel.nugget)==0) {
      if(linear.ID.coords) {
         out[object$p+3] <- out[object$p+1]*out[object$p+3]
         names(out)[object$p+3] <- "tau^2"
         out[object$p+4] <- out[object$p+1]*out[object$p+4]
         names(out)[object$p+4] <- "omega^2"
      } else {
         out[object$p+3] <- out[object$p+1]*out[object$p+3]
         names(out)[object$p+3] <- "tau^2"
      }
   }
   return(out)
}

##' @title Plot of a predicted surface
##' @description \code{plot.pred.PrevMap} displays predictions obtained from \code{\link{spatial.pred.linear.MLE}}, \code{\link{spatial.pred.linear.Bayes}},\code{\link{spatial.pred.binomial.MCML}}, \code{\link{spatial.pred.binomial.Bayes}} and \code{\link{spatial.pred.poisson.MCML}}.
##' @param x an object of class "PrevMap".
##' @param type a character indicating the type of prediction to display: 'prevalence','odds', 'logit' or 'probit' for binomial models; "log" or "exponential" for Poisson models. Default is \code{NULL}.
##' @param summary character indicating which summary to display: 'predictions','quantiles', 'standard.errors' or 'exceedance.prob'; default is 'predictions'. If \code{summary="exceedance.prob"}, the argument \code{type} is ignored.
##' @param ... further arguments passed to \code{\link{plot}} of the 'raster' package.
##' @method plot pred.PrevMap
##' @importFrom raster rasterFromXYZ
##' @importFrom methods getMethod signature
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
plot.pred.PrevMap <- function(x,type=NULL,summary="predictions",...) {
	if(class(x)!="pred.PrevMap") stop("x must be of class pred.PrevMap")
	if(length(type)>0 && 
	   any(type==c("prevalence","odds","logit","probit","log","exponential"))==FALSE) {
	   stop("type must be 'prevalence','odds', 'logit' or 'probit' for the binomial modl, and 'log' or 'exponential' for the Poisson model.")
	}
	
	if(length(type)>0 & summary=="exceedance.prob") {
		   warning("the argument 'type' is ignored when summary='exceedance.prob'")
    }
	
	if(length(type)==0 && summary!="exceedance.prob") {
	   stop("type must be specified")	
	}
	
	if(summary !="exceedance.prob") {
   	     
   	   if(any(summary==c("predictions",
   	     "quantiles","standard.errors","exceedance.prob"))==FALSE) {
		   stop("summary must be 'predictions','quantiles',
		        'standard.errors' or 'exceedance.prob'")
   	   }
   	   
   	   r <- rasterFromXYZ(cbind(x$grid,x[[type]][[summary]]))
	} else {
	   r <- rasterFromXYZ(cbind(x$grid,x[[summary]]))       
	}
	getMethod('plot',signature=signature(x='Raster', y='ANY'))(r,...)
}

##' @title Contour plot of a predicted surface 
##' @description \code{plot.pred.PrevMap} displays contours of predictions obtained from \code{\link{spatial.pred.linear.MLE}}, \code{\link{spatial.pred.linear.Bayes}},\code{\link{spatial.pred.binomial.MCML}} and \code{\link{spatial.pred.binomial.Bayes}}.
##' @param x an object of class "pred.PrevMap".
##' @param type a character indicating the type of prediction to display: 'prevalence', 'odds', 'logit' or 'probit'.
##' @param ... further arguments passed to \code{\link{contour}}.
##' @param summary character indicating which summary to display: 'predictions','quantiles', 'standard.errors' or 'exceedance.prob'; default is 'predictions'. If \code{summary="exceedance.prob"}, the argument \code{type} is ignored.
##' @method contour pred.PrevMap
##' @importFrom raster rasterFromXYZ
##' @importFrom methods getMethod signature
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
contour.pred.PrevMap <- function(x,type=NULL,summary="predictions",...) {
	
   if(class(x)!="pred.PrevMap") stop("x must be of class pred.PrevMap")
	
   if(length(type)>0 && 
	   any(type==c("prevalence","odds","logit","probit"))==FALSE) {
	   stop("type must be 'prevalence','odds', 'logit' or 'probit'")
   }
	
   if(length(type)>0 & summary=="exceedance.prob") {
		warning("the argument 'type' is ignored when 
		         summary='exceedance.prob'")
   }
	
   if(length(type)==0 && summary!="exceedance.prob") {
	   stop("type must be specified")	
   }
	
   if(summary !="exceedance.prob") {
      if(any(summary==c("predictions","quantiles",
   	        	            "standard.errors"))==FALSE) {
		   stop("summary must be 'predictions','quantiles', 
		         'standard.errors' or 'exceedance.prob'")
   	  }
   	  r <- rasterFromXYZ(cbind(x$grid,x[[type]][[summary]]))   	   
   } else {
	  r <- rasterFromXYZ(cbind(x$grid,x[[summary]]))
   }
   getMethod('contour',signature=signature(x='RasterLayer'))(r,...)
}

##' @title Summarizing Bayesian model fits
##' @description \code{summary} method for the class "Bayes.PrevMap" that computes the posterior mean, median, mode and high posterior density intervals using samples from Bayesian fits.
##' @param object an object of class "Bayes.PrevMap" obatained as result of a call to \code{\link{binomial.logistic.Bayes}} or \code{\link{linear.model.Bayes}}.
##' @param hpd.coverage value of the coverage of the high posterior density intervals; default is \code{0.95}.
##' @param ... further arguments passed to or from other methods.
##' @return A list with the following values
##' @return \code{linear}: logical value that is \code{TRUE} if a linear model was fitted and \code{FALSE} otherwise.
##' @return \code{binary}: logical value that is \code{TRUE} if a binary model was fitted and \code{FALSE} otherwise.
##' @return \code{probit}: logical value that is \code{TRUE} if a binary model with probit link function was fitted and \code{FALSE} if with logistic link function.
##' @return \code{ck}: logical value that is \code{TRUE} if a low-rank approximation was fitted and \code{FALSE} otherwise.
##' @return \code{beta}: matrix of the posterior summaries for the regression coefficients.
##' @return \code{sigma2}: vector of the posterior summaries for \code{sigma2}.
##' @return \code{phi}: vector of the posterior summaries for \code{phi}.
##' @return \code{tau2}: vector of the posterior summaries for \code{tau2}.
##' @return \code{call}: matched call.
##' @return \code{kappa}: fixed value of the shape paramter of the Matern covariance function.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method summary Bayes.PrevMap
##' @export
summary.Bayes.PrevMap <- function(object,hpd.coverage=0.95,...) {
	hpd <- function(x, conf){
	   conf <- min(conf, 1-conf)
	   n <- length(x)
	   nn <- round( n*conf )
	   x <- sort(x)
	   xx <- x[ (n-nn+1):n ] - x[1:nn]
   	   m <- min(xx)
	   nnn <- which(xx==m)[1]
	   return( c( x[ nnn ], x[ n-nn+nnn ] ) )
    }
    post.mode <- function(x) {   
    	   d <- density(x) 	
    	   d$x[which.max(d$y)]
    }
	if(class(object)!="Bayes.PrevMap") stop("object must be of class 'PrevMap'.")
	res <- list()
	
	check.binary <- all(object$y==0 | object$y==1)

	res$linear <- FALSE
	res$binary <- FALSE
	res$probit <- FALSE
	
	if (check.binary){		
		res$binary <- TRUE
		if(substr(object$call[1],8,13)=="probit") {
		   res$probit <- TRUE
		}
    } else if (length(object$units.m)==0 & !check.binary) {
		res$linear <- TRUE
	}
	
	
	if(length(object$knots)>0) {
		res$ck <- TRUE
	} else {
		res$ck <- FALSE
	}
	p <- ncol(object$D)
	res$beta <- matrix(NA,nrow=p,ncol=6)	
	for(i in 1:p) {	
		res$beta[i,] <- c(mean(as.matrix(object$estimate[,1:p])[,i]),
		                          median(as.matrix(object$estimate[,1:p])[,i]),
		                          post.mode(as.matrix(object$estimate[,1:p])[,i]),
	                              sd(as.matrix(object$estimate[,1:p])[,i]),
	                              hpd(as.matrix(object$estimate[,1:p])[,i],hpd.coverage))
	}
	names.summaries <- c("Mean","Median","Mode","StdErr", paste("HPD",(1-hpd.coverage)/2),
	                                      paste("HPD",1-(1-hpd.coverage)/2))
	rownames(res$beta) <- colnames(object$D)
	colnames(res$beta) <- names.summaries
    res$sigma2 <- 	  c(mean(object$estimate[,"sigma^2"]),median(object$estimate[,"sigma^2"]),
                                 post.mode(object$estimate[,"sigma^2"]),
	                             sd(object$estimate[,"sigma^2"]),
	                             hpd(object$estimate[,"sigma^2"],hpd.coverage))   
	names(res$sigma2) <- names.summaries                                        
    res$phi <- 	  c(mean(object$estimate[,"phi"]),median(object$estimate[,"phi"]),
                                 post.mode(object$estimate[,"phi"]),
	                             sd(object$estimate[,"phi"]),
	                             hpd(object$estimate[,"phi"],hpd.coverage))	
    names(res$phi) <- names.summaries     	                             
    if(ncol(object$estimate)==p+3) {
    	res$tau2 <- 	  c(mean(object$estimate[,"tau^2"]),
    	                         median(object$estimate[,"tau^2"]),
                             post.mode(object$estimate[,"tau^2"]),
	                             sd(object$estimate[,"tau^2"]),
	                             hpd(object$estimate[,"tau^2"],hpd.coverage))
	    names(res$tau2) <- names.summaries                              
    }	       
    res$call <- object$call 
    res$kappa <- object$kappa        
    class(res) <- "summary.Bayes.PrevMap"
    return(res)                                                      
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method print summary.Bayes.PrevMap
##' @export
print.summary.Bayes.PrevMap <- function(x,...) {
	if(x$linear) {
       cat("Bayesian geostatistical linear model \n")
    } else if(x$probit) {
       cat("Bayesian binary geostatistical probit model \n") 
    } else {	
       cat("Bayesian binomial geostatistical logistic model \n")    	
    }
	if(x$ck) {
	  cat("(low-rank approximation) \n")
	} 	
	cat("Call: \n")
	print(x$call)
	cat("\n")
	print(x$beta)
	cat("\n")
	if(x$ck) {
	   cat("Matern kernel parameters (kappa=",
	      x$kappa,") \n \n",sep="")
	} else{
       cat("Covariance parameters Matern function (kappa=",
	      x$kappa,") \n \n",sep="")
	}      
	if(length(x$tau2) > 0) {
	   tab <- rbind(x$sigma2,x$phi,x$tau2)
	   rownames(tab) <- c("sigma^2","phi","tau^2")
	   print(tab) 
	} else {
	   tab <- rbind(x$sigma2,x$phi)
	   rownames(tab) <- c("sigma^2","phi")
	   print(tab) 
	}	   
	cat("\n")
    cat("Legend: \n")
    if(x$ck) {
        cat("sigma^2 = variance of the iid zero-mean Gaussian variables \n")	
    } else {
    	cat("sigma^2 = variance of the Gaussian process \n")
    }
    cat("phi = scale of the spatial correlation \n")
    if(length(x$tau2) > 0) cat("tau^2 = variance of the nugget effect \n")
}

##' @title Plot of the autocorrelgram for posterior samples 
##' @description Plots the autocorrelogram for the posterior samples of the model parameters and spatial random effects.
##' @param object an object of class 'Bayes.PrevMap'.
##' @param param a character indicating for which component of the model the autocorrelation plot is required: \code{param="beta"} for the regression coefficients; \code{param="sigma2"} for the variance of the spatial random effect; \code{param="phi"} for the scale parameter of the Matern correlation function; \code{param="tau2"} for the variance of the nugget effect; \code{param="S"} for the spatial random effect.
##' @param component.beta if \code{param="beta"}, \code{component.beta} is a numeric value indicating the component of the regression coefficients; default is \code{NULL}.
##' @param component.S if \code{param="S"}, \code{component.S} can be a numeric value indicating the component of the spatial random effect, or set equal to \code{"all"} if the autocorrelgram should be plotted for all the components. Default is \code{NULL}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
autocor.plot <- function(object,param,component.beta=NULL,component.S=NULL) {
    p <- ncol(object$D)
	if(class(object)!="Bayes.PrevMap") stop("object must be of class 'PrevMap'.")
	if(any(param==c("beta","sigma2","phi","tau2","S"))==FALSE) stop("param must be equal to 'beta', 'sigma2','phi' or 'S'.")		
	if(param=="beta" && length(component.beta)==0) stop("if param='beta', component.beta must be provided.")
	if(param=="beta" && component.beta > p) stop("wrong value for component.beta")
	if(param=="beta" && component.beta <= 0) stop("component.beta must be positive")
	if(param=="tau2" && ncol(object$estimate)!=p+3) stop("the nugget effect was not included in the model.")	
	if(param=="S" && length(component.S)==0) stop("if param='S', component.S must be provided.")	
	if(is.character(component.S) && component.S !="all") stop("component.S must be positive, or equal to 'all' if the autocorrelogram plot is required for all the components of the spatial random effect.")
	if(is.numeric(component.S) && component.S <= 0) stop("component.S must be positive, or equal to 'all' if the autocorrelogram plot is required for all the components of the spatial random effect.")
	if(length(component.S) > 0 && is.character(component.S)==FALSE && is.numeric(component.S)==FALSE) stop("component.S must be positive, or equal to 'all' if the autocorrelogram plot is required for all the components of the spatial random effect.")
    if(param=="S" && is.character(component.S) && component.S=="all")  {
    	   acf.plot <- acf(object$S[,1],plot=FALSE)
       plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
                 ylim=c(-1,1))                
       for(i in 2:ncol(object$S)) {
          acf.plot <- acf(object$S[,i],plot=FALSE)
          lines(acf.plot$lag,acf.plot$acf)
       }
    	   abline(h=0,lty="dashed",col=2)
    } else if(param=="S" && is.numeric(component.S)) {
    	   acf(object$S[,component.S],main=paste("Spatial random effect: comp. n.",component.S,sep=""))
    } else if(param=="beta") {
    	   acf(as.matrix(object$estimate[,1:p])[,component.beta],main=paste(colnames(object$D)[component.beta]))
    } else if(param=="sigma2") {
    	   acf(object$estimate[,"sigma^2"],main=expression(sigma^2))
    } else if(param=="phi") {
       acf(object$estimate[,"phi"],main=expression(phi))
    } else if(param=="tau2") {
       acf(object$estimate[,"tau^2"],main=expression(tau^2))
    }
}

##' @title Trace-plots for posterior samples 
##' @description Displays the trace-plots for the posterior samples of the model parameters and spatial random effects.
##' @param object an object of class 'Bayes.PrevMap'.
##' @param param a character indicating for which component of the model the density plot is required: \code{param="beta"} for the regression coefficients; \code{param="sigma2"} for the variance of the spatial random effect; \code{param="phi"} for the scale parameter of the Matern correlation function; \code{param="tau2"} for the variance of the nugget effect; \code{param="S"} for the spatial random effect.
##' @param component.beta if \code{param="beta"}, \code{component.beta} is a numeric value indicating the component of the regression coefficients; default is \code{NULL}.
##' @param component.S if \code{param="S"}, \code{component.S} can be a numeric value indicating the component of the spatial random effect. Default is \code{NULL}.
##' @param ... additional parameters to pass to \code{\link{density}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
trace.plot <- function(object,param,component.beta=NULL,component.S=NULL) {
	    p <- ncol(object$D)
	if(class(object)!="Bayes.PrevMap") stop("object must be of class 'PrevMap'.")
	if(any(param==c("beta","sigma2","phi","tau2","S"))==FALSE) stop("param must be equal to 'beta', 'sigma2','phi' or 'S'.")		
	if(param=="beta" && length(component.beta)==0) stop("if param='beta', component.beta must be provided.")
	if(param=="beta" && component.beta > p) stop("wrong value for component.beta")
	if(param=="beta" && component.beta <= 0) stop("component.beta must be positive")
	if(param=="tau2" && ncol(object$estimate)!=p+3) stop("the nugget effect was not included in the model")
    if(param=="S" && length(component.S)==0) stop("if param='S', component.S must be provided.")	
	if(is.numeric(component.S) && component.S <= 0) stop("component.S must be positive integer.")
    if(param=="S") {
       plot(object$S[,component.S],type="l",xlab="Iteration",
    	   main=paste("Spatial random effect: comp. n.",component.S,sep=""),
    	   ylab="")
    } else if(param=="beta") {
       plot(as.matrix(object$estimate[,1:p])[,component.beta],main=paste(colnames(object$D)[component.beta]),
    	   type="l",xlab="Iteration",ylab="")
    } else if(param=="sigma2") {
       plot(object$estimate[,"sigma^2"],main=expression(sigma^2),type="l",xlab="Iteration",ylab="")
    } else if(param=="phi") {
       plot(object$estimate[,"phi"],main=expression(phi),xlab="Iteration",type="l",ylab="")
    } else if(param=="tau2") {
       plot(object$estimate[,"tau^2"],main=expression(tau^2),xlab="Iteration",type="l",ylab="")
    }
}

##' @title Density plot for posterior samples
##' @description Plots the autocorrelogram for the posterior samples of the model parameters and spatial random effects.
##' @param object an object of class 'Bayes.PrevMap'.
##' @param param a character indicating for which component of the model the density plot is required: \code{param="beta"} for the regression coefficients; \code{param="sigma2"} for the variance of the spatial random effect; \code{param="phi"} for the scale parameter of the Matern correlation function; \code{param="tau2"} for the variance of the nugget effect; \code{param="S"} for the spatial random effect.
##' @param component.beta if \code{param="beta"}, \code{component.beta} is a numeric value indicating the component of the regression coefficients; default is \code{NULL}.
##' @param component.S if \code{param="S"}, \code{component.S} can be a numeric value indicating the component of the spatial random effect. Default is \code{NULL}.
##' @param hist logical; if \code{TRUE} a histrogram is added to density plot.
##' @param ... additional parameters to pass to \code{\link{density}}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
dens.plot <- function(object,param,component.beta=NULL,
                                          component.S=NULL,hist=TRUE,...) {
    p <- ncol(object$D)
	if(class(object)!="Bayes.PrevMap") stop("object must be of class 'Bayes.PrevMap'.")
	if(any(param==c("beta","sigma2","phi","tau2","S"))==FALSE) stop("param must be equal to 'beta', 'sigma2','phi' or 'S'.")		
	if(param=="beta" && length(component.beta)==0) stop("if param='beta', component.beta must be provided.")
	if(param=="beta" && component.beta > p) stop("wrong value for component.beta")
	if(param=="beta" && component.beta <= 0) stop("component.beta must be positive")
	if(param=="tau2" && ncol(object$estimate)!=p+3) stop("the nugget effect was not included in the model.")	
	if(param=="S" && length(component.S)==0) stop("if param='S', component.S must be provided.")	
	if(is.character(component.S) && component.S !="all") stop("component.S must be a positive integer.")
	if(is.numeric(component.S) && component.S <= 0) stop("component.S must be positive integer.")
    if(param=="S") {
    	    if(hist) {
        	    dens <- density(object$S[,component.S],...)
            hist(object$S[,component.S],
                  main=paste("Spatial random effect: comp. n.",component.S,sep=""),
                  prob=TRUE,xlab="",ylim=c(0,max(dens$y)))          
            lines(dens,lwd=2)
        } else {
        	    dens <- density(object$S[,component.S],...)
            plot(dens,
            main=paste("Spatial random effect: comp. n.",component.S,sep=""),
            lwd=2,xlab="")
        }  
    } else if(param=="beta") {
        if(hist) {
        	   dens <- density(as.matrix(object$estimate[,1:p])[,component.beta],...)
           hist(as.matrix(object$estimate[,1:p])[,component.beta],
           main=paste(colnames(object$D)[component.beta]),
           prob=TRUE,xlab="",ylim=c(0,max(dens$y)))          
           lines(dens,lwd=2)
        } else {
        	   dens <- density(as.matrix(object$estimate[,1:p])[,component.beta],...)
           plot(dens,main=paste(colnames(object$D)[component.beta]),lwd=2,xlab="")
        }   
    } else if(param=="sigma2") {
        if(hist) {
           dens <- density(object$estimate[,"sigma^2"],...)
           hist(object$estimate[,"sigma^2"],main=expression(sigma^2),prob=TRUE,
           xlab="",ylim=c(0,max(dens$y)))          
           lines(dens,lwd=2)
        } else {
           dens <- density(object$estimate[,"sigma^2"],...)
           plot(dens,main=expression(sigma^2),lwd=2,xlab="")
        }     
    } else if(param=="phi") {
        if(hist) {
           dens <- density(object$estimate[,"phi"],...)
           hist(object$estimate[,"phi"],main=expression(phi),prob=TRUE,
           xlab="",ylim=c(0,max(dens$y)))          
           lines(dens,lwd=2)
        } else {
           dens <- density(object$estimate[,"phi"],...)
           plot(dens,main=expression(phi),lwd=2,xlab="")
        }  
    } else if(param=="tau2") {
        if(hist) {
           dens <- density(object$estimate[,"tau^2"],...)
           hist(object$estimate[,"tau^2"],main=expression(tau^2),prob=TRUE,
           xlab="",ylim=c(0,max(dens$y)))          
           lines(dens,lwd=2)
        } else {
           dens <- density(object$estimate[,"tau^2"],...)
           plot(dens,main=expression(tau^2),lwd=2,xlab="")
        }  
    }
}

##' @title Profile likelihood for the shape parameter of the Matern covariance function
##' @description This function plots the profile likelihood for the shape parameter of the Matern covariance function used in the linear Gaussian model. It also computes confidence intervals of coverage \code{coverage} by interpolating the profile likelihood with a spline and using the asymptotic distribution of a chi-squared with one degree of freedom.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model.
##' @param set.kappa a vector indicating the set values for evluation of the profile likelihood.
##' @param fixed.rel.nugget a value for the relative variance \code{nu2} of the nugget effect, that is then treated as fixed. Default is \code{NULL}.
##' @param start.par starting values for the scale parameter \code{phi} and the relative variance of the nugget effect \code{nu2}; if \code{fixed.rel.nugget} is provided, then a starting value for \code{phi} only should be provided.
##' @param coverage a value between 0 and 1 indicating the coverage of the confidence interval based on the interpolated profile liklelihood for the shape parameter. Default is \code{coverage=NULL} and no confidence interval is then computed.
##' @param plot.profile logical; if \code{TRUE} the computed profile-likelihood is plotted together with the interpolating spline. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return The function returns an object of class 'shape.matern' that is a list with the following components
##' @return \code{set.kappa} set of values of the shape parameter used to evaluate the profile-likelihood.
##' @return \code{val.kappa} values of the profile likelihood. 
##' @return If a value for \code{coverage} is specified, the list also contains \code{lower}, \code{upper} and \code{kappa.hat} that correspond to the lower and upper limits of the confidence interval, and the maximum likelihood estimate for the shape parameter, respectively.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom geoR varcov.spatial
##' @export
shape.matern <- function(formula,coords,data,set.kappa,fixed.rel.nugget=NULL,start.par,
              coverage=NULL,plot.profile=TRUE,messages=TRUE) {              	
	start.par <- as.numeric(start.par)
	if(length(coverage)>0 && (coverage <= 0 | coverage >=1)) stop("'coverage' must be between 0 and 1.")
	mf <- model.frame(formula,data=data)
	if(any(is.na(data))) stop("missing data are not accepted.")	
	y <- as.numeric(model.response(mf))
	n <- length(y)
	D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
	coords <- as.matrix(model.frame(coords,data))
	if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("wrong set of coordinates.")
	U <- dist(coords)
	p <- ncol(D)
	if((length(start.par)!= 2 & length(fixed.rel.nugget)==0) |
	   (length(start.par)!= 1 & length(fixed.rel.nugget)>0)) stop("wrong length of start.cov.pars")
	if(any(set.kappa < 0)) stop("if log.kappa=FALSE, then set.kappa must have positive components.")  
	log.lik.profile <- function(phi,kappa,nu2) {
	      V <- varcov.spatial(dists.lowertri=U,cov.pars=c(1,phi),kappa=kappa,
	              nugget=nu2)$varcov
	      V.inv <- solve(V)
	      beta.hat <- as.numeric(solve(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%y)
	      diff.beta.hat <- y-D%*%beta.hat
	      sigma2.hat <- as.numeric(t(diff.beta.hat)%*%V.inv%*%diff.beta.hat/n)
	      ldet.V <- determinant(V)$modulus
	      as.numeric(-0.5*(n*log(sigma2.hat)+ldet.V))
    }
    start.par <- log(start.par)
    val.kappa <- rep(NA,length(set.kappa))
    if(length(fixed.rel.nugget) == 0) {
       for(i in 1:length(set.kappa)) {
          val.kappa[i] <- -nlminb(start.par,function(x) 
       	                         -log.lik.profile(exp(x[1]),set.kappa[i],exp(x[2])))$objective
       	  if(messages) cat("Evalutation for kappa=",set.kappa[i],"\r",sep="")                         
          flush.console()
       }    	
    } else {    	
        for(i in 1:length(set.kappa)) {      	         
          val.kappa[i] <- -nlminb(start.par,function(x) 
       	                         -log.lik.profile(exp(x),set.kappa[i],fixed.rel.nugget))$objective
       	  if(messages) cat("Evalutation for kappa=",set.kappa[i],"\r",sep="")                         
          flush.console()
       }  	
    }
    out <- list()
    out$set.kappa <- set.kappa
    out$val.kappa <- val.kappa
    if(length(coverage) > 0) {
      	f <- splinefun(set.kappa,val.kappa)
    	    estim.kappa <- optimize(function(x)-f(x),
                                lower=min(set.kappa),upper=max(set.kappa))
        max.log.lik <- -estim.kappa$objective                       
        par.set <- seq(min(set.kappa),max(set.kappa),length=20)
        k <- qchisq(coverage,df=1)
        par.val <- 2*(max.log.lik-f(par.set))-k
        done <- FALSE
        if(par.val[1] < 0 & par.val[length(par.val)] > 0) {
        	   stop("smaller value of the shape parameter should be included when computing the profile likelihood")
        } else if(par.val[1] > 0 & par.val[length(par.val)] < 0) {
        	   stop("larger value of the shape parameter should be included when computing the profile likelihood")
        } else if(par.val[1] < 0 & par.val[length(par.val)] < 0) {
           stop("both smaller and larger value of the parameter should be included when computing the profile likelihood")        	
        }
        i <- 1      
        lower.int <- c(NA,NA)
        upper.int <- c(NA,NA)        
        while(!done) {
        	    i <- i+1           	
           	if(sign(par.val[i-1]) != sign(par.val[i]) & sign(par.val[i-1])==1) {
           		lower.int[1] <- par.set[i-1]
           		lower.int[2] <- par.set[i+1]
           	}           	
           	if(sign(par.val[i-1]) != sign(par.val[i]) & sign(par.val[i-1])==-1) {
           		upper.int[1] <- par.set[i-1]
           		upper.int[2] <- par.set[i+1]
           		done <- TRUE
           	}       
        }
        out$lower <- uniroot(function(x) 2*(max.log.lik-f(x))-k,
        lower=lower.int[1],upper=lower.int[2])$root                        
        out$upper <- uniroot(function(x) 2*(max.log.lik-f(x))-k,
        lower=upper.int[1],upper=upper.int[2])$root 
        out$kappa.hat <- estim.kappa$minimum
      } 
      if(plot.profile) {
    	  plot(set.kappa,val.kappa,type="l",xlab=expression(kappa),ylab="log-likelihood",
                    main=expression("Profile likelihood for"~kappa))
          if(length(coverage) > 0) {          
             a <- seq(min(par.set),max(par.set),length=1000)
             b <- sapply(a, function(x) f(x))
             lines(a,b,type="l",col=2)
             abline(h=-(k/2-max.log.lik),lty="dashed")
             abline(v=out$lower,lty="dashed")
             abline(v=out$upper,lty="dashed")
       }
    }          
    class(out) <- "shape.matern"
    return(out)
}

##' @title Plot of the profile likelihood for the shape parameter of the Matern covariance function
##' @description This function plots the profile likelihood for the shape parameter of the Matern covariance function using the output from \code{\link{shape.matern}} function.
##' @param x an object of class 'shape.matern' obtained as result of a call to \code{\link{shape.matern}}
##' @param plot.spline logical; if \code{TRUE} an interpolating spline of the profile likelihood is added to the plot. 
##' @param ... further arguments passed to \code{\link{plot}}.
##' @seealso \code{\link{shape.matern}}
##' @return The function does not return any value.
##' @method plot shape.matern
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
plot.shape.matern <- function(x,plot.spline=TRUE,...) {
	if(class(x)!="shape.matern")  stop("x must be of class 'shape.matern'")
	
	plot.list <- list(...)
	plot.list$x <- x$set.kappa
	plot.list$y <- x$val.kappa
	if(length(plot.list$main)==0) plot.list$main <- expression("Profile likelihood for"~kappa)
	if(length(plot.list$xlab)==0) plot.list$xlab <- expression(kappa)
	if(length(plot.list$ylab)==0) plot.list$ylab <- "log-likelihood"
	if(length(plot.list$type)==0) plot.list$type <- "l"
	do.call(plot,plot.list)
	if(plot.spline) {
	   f <- splinefun(x$set.kappa,x$val.kappa)                     
       par.set <- seq(min(x$set.kappa),max(x$set.kappa),length=100)
       par.val <- f(par.set)
       lines(par.set,par.val,col=2)
	}  	  
}

##' @title Adjustment factor for the variance of the convolution of Gaussian noise
##' @description This function computes the multiplicative constant used to adjust the value of \code{sigma2} in the low-rank approximation of a Gaussian process. 
##' @param knots.dist a matrix of the distances between the observed coordinates and the spatial knots.
##' @param phi scale parameter of the Matern covariance function.
##' @param kappa shape parameter of the Matern covariance function.
##' @details Let \eqn{U} denote the \eqn{n} by \eqn{m} matrix of the distances between the \eqn{n} observed coordinates and \eqn{m} pre-defined spatial knots. This function computes the following quantity 
##' \deqn{\frac{1}{n}\sum_{i=1}^n \sum_{j=1}^m K(u_{ij}; \phi, \kappa)^2,}
##' where \eqn{K(.; \phi, \kappa)} is the Matern kernel (see \code{\link{matern.kernel}}) and \eqn{u_{ij}} is the distance between the \eqn{i}-th sampled location and the \eqn{j}-th spatial knot. 
##' @return A value corresponding to the adjustment factor for \code{sigma2}.
##' @seealso \code{\link{matern.kernel}}, \code{\link{pdist}}. 
##' @examples
##' set.seed(1234)
##' # Observed coordinates
##' n <- 100
##' coords <- cbind(runif(n),runif(n))
##' 
##' # Spatial knots
##' knots <- expand.grid(seq(-0.2,1.2,length=5),seq(-0.2,1.2,length=5))
##' 
##' # Distance matrix
##' knots.dist <- as.matrix(pdist(coords,knots))
##' 
##' adjust.sigma2(knots.dist,0.1,2)
##' 
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export

adjust.sigma2 <- function(knots.dist,phi,kappa) {
	K <- matern.kernel(knots.dist,2*sqrt(kappa)*phi,kappa)
	out <- mean(apply(K,1,function(r) sum(r^2)))
	return(out)
}

##' @title Spatially discrete sampling
##' @description Draws a sub-sample from a set of units spatially located irregularly over some defined geographical region by imposing a minimum distance between any two sampled units.  
##' @param xy.all set of locations from which the sample will be drawn.
##' @param n size of required sample.
##' @param delta minimum distance between any two locations in preliminary sample.
##' @param k number of locations in preliminary sample to be replaced by nearest neighbours of other preliminary sample locations in final sample (must be between 0 and \code{n/2}).
##'
##' @details  To draw a sample of size \code{n}  from a population of spatial locations \eqn{X_{i}  : i  = 1,\ldots,N}, with the property that the distance between any two sampled locations is at least \code{delta}, the function implements the following algorithm.
##' \itemize{
##' \item{Step 1.} Draw an initial sample of size \code{n}  completely at random and call this \eqn{x_{i}  : i  = 1,\dots, n}.
##' \item{Step 2.} Set \eqn{i  = 1} and calculate the minimum, \eqn{d_{\min}}, of the distances from \eqn{x_{i}}  to all other \eqn{x_{j}}  in the initial sample.
##' \item{Step 3.} If \eqn{d_{\min} \ge \delta}, increase \eqn{i}  by 1 and return to step 2 if \eqn{i \le n}, otherwise stop.
##' \item{Step 4.} If \eqn{d_{\min} < \delta}, draw an integer \eqn{j}  at random from \eqn{1,  2,\ldots,N}, set \eqn{x_{i}  = X_{j}}  and return to step 3.
##' }
##' Samples generated in this way will exhibit a more regular spatial arrangement than would a random sample of the same size. The degree of regularity achievable will be influenced by the spatial arrangement of the population \eqn{X_{i}  : i  = 1,\ldots,N}, the specified value of \code{delta}  and the sample size \code{n}. For any given population, if \code{n}  and/or \code{delta}  are too large, a sample of the required size with the distance between any two sampled locations at least \code{delta} will not be achievable; the suggested solution is then to run the algorithm with a smaller value of \code{delta}.
##'
##' \bold{Sampling close pairs of points}.
##'  For some purposes, it is desirable that a spatial sampling scheme include pairs of closely spaced points. In this case, the above algorithm requires the following additional steps to be taken.
##' Let \code{k}  be the required number of close pairs.
##' \itemize{
##' \item{Step 5.} Set \eqn{j  = 1} and draw a random sample of size 2 from the integers \eqn{1,  2,\ldots,n}, say \eqn{(i_{1}, i_{2} )}.
##' \item{Step 6.} Find the integer \eqn{r}  such that the distances from \eqn{x_{i_{1}}}  to \eqn{X_{r}} is the minimum of all \eqn{N-1} distances from \eqn{x_{i_{1}}}  to the \eqn{X_{j}}.
##' \item{Step 7.}  Replace \eqn{x_{i_{2}}}  by \eqn{X_{r}}, increase \eqn{i}  by 1 and return to step 5 if \eqn{i \le k}, otherwise stop.
##' }
##' 
##' @return A matrix of dimension \code{n} by 2 containing the final sampled locations.
##'
##' @examples
##' x<-0.015+0.03*(1:33)
##' xall<-rep(x,33)
##' yall<-c(t(matrix(xall,33,33)))
##' xy<-cbind(xall,yall)+matrix(-0.0075+0.015*runif(33*33*2),33*33,2)
##' par(pty="s",mfrow=c(1,2))
##' plot(xy[,1],xy[,2],pch=19,cex=0.25,xlab="Easting",ylab="Northing",
##'    cex.lab=1,cex.axis=1,cex.main=1)
##' 
##' set.seed(15892) 
##' # Generate spatially random sample
##' xy.sample<-xy[sample(1:dim(xy)[1],50,replace=FALSE),] 
##' points(xy.sample[,1],xy.sample[,2],pch=19,col="red")
##' points(xy[,1],xy[,2],pch=19,cex=0.25)
##' plot(xy[,1],xy[,2],pch=19,cex=0.25,xlab="Easting",ylab="Northing",
##'    cex.lab=1,cex.axis=1,cex.main=1)
##'
##' set.seed(15892) 
##' # Generate spatially regular sample
##' xy.sample<-discrete.sample(xy,50,0.08) 
##' points(xy.sample[,1],xy.sample[,2],pch=19,col="red")
##' points(xy[,1],xy[,2],pch=19,cex=0.25)
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##' @export

discrete.sample<-function(xy.all,n,delta,k=0) {
   deltasq<-delta*delta
   N<-dim(xy.all)[1]
   index<-1:N
   index.sample<-sample(index,n,replace=FALSE)
   xy.sample<-xy.all[index.sample,]
   for (i in 2:n) {
      dmin<-0
      while (dmin<deltasq) {
         take<-sample(index,1)
         dvec<-(xy.all[take,1]-xy.sample[,1])^2+(xy.all[take,2]-xy.sample[,2])^2
         dmin<-min(dvec)
         }
      xy.sample[i,]<-xy.all[take,]
      }
   if (k>0) {
   	  if(k > n/2) stop("k must be between 0 and n/2.")
      take<-matrix(sample(1:n,2*k,replace=FALSE),k,2)
      for (j in 1:k) {
         take1<-take[j,1]; take2<-take[j,2]
         xy1<-c(xy.sample[take1,])
         dvec<-(xy.all[,1]-xy1[1])^2+(xy.all[,2]-xy1[2])^2
         neighbour<-order(dvec)[2]
         xy.sample[take2,]<-xy.all[neighbour,]
         }
      }
   return(xy.sample)
} 

##' @title Spatially continuous sampling
##' @description Draws a sample of spatial locations within a spatially continuous polygonal sampling region.
##' @param poly boundary of a polygon.
##' @param n number of events.
##' @param delta minimum permissible distance between any two events in preliminary sample.
##' @param k number of locations in preliminary sample to be replaced by near neighbours of other preliminary sample locations in final sample (must be between 0 and \code{n/2})
##' @param rho maximum distance between close pairs of locations in final sample.  
##'
##' @return A matrix of dimension \code{n} by 2 containing event locations.
##'  
##' @details  To draw a sample of size \code{n}  from a spatially continuous region \eqn{A}, with the property that the distance between any two sampled locations is at least \code{delta}, the following algorithm is used.
##' \itemize{
##' \item{Step 1.} Set \eqn{i  = 1} and generate a point \eqn{x_{1}}  uniformly distributed on \eqn{A}.
##' \item{Step 2.} Increase \eqn{i}  by 1, generate a point \eqn{x_{i}}  uniformly distributed on \eqn{A} and calculate the minimum, \eqn{d_{\min}}, of the distances from \eqn{x_{i}} to all \eqn{x_{j}: j < i }.
##' \item{Step 3.} If \eqn{d_{\min} \ge \delta}, increase \eqn{i}  by 1 and return to step 2 if \eqn{i \le n}, otherwise stop;
##' \item{Step 4.} If \eqn{d_{\min} < \delta}, return to step 2 without increasing \eqn{i}.
##' }
##'
##' \bold{ Sampling close pairs of points.}  For some purposes, it is desirable that a spatial sampling scheme include pairs of closely spaced points. In this case, the above algorithm requires the following additional steps to be taken.
##'  Let \code{k}  be the required number of close pairs. Choose a value \code{rho}  such that a close pair  of points will be a pair of points separated by a distance of at most \code{rho}.
##' \itemize{
##' \item{Step 5.} Set \eqn{j  = 1} and draw a random sample of size 2 from the integers \eqn{1,2,\ldots,n}, say \eqn{(i_{1}; i_{2})};
##' \item{Step 6.} Replace \eqn{x_{i_{1}}} by \eqn{x_{i_{2}} + u} , where \eqn{u}  is uniformly distributed on the disc with centre \eqn{x_{i_{2}}} and radius \code{rho}, increase \eqn{i} by 1 and return to step 5 if \eqn{i \le k}, otherwise stop.	
##' }	
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##' @examples
##' library(geoR)
##' data(parana)
##' poly<-parana$borders
##' poly<-matrix(c(poly[,1],poly[,2]),dim(poly)[1],2,byrow=FALSE)
##' set.seed(5871121)
##' 
##' # Generate spatially regular sample
##' xy.sample<-continuous.sample(poly,100,30) 
##' plot(poly,type="l",xlab="X",ylab="Y")
##' points(xy.sample,pch=19,cex=0.5)
##'
##' @importFrom splancs csr
##' @export
continuous.sample<-function(poly,n,delta,k=0,rho=NULL) {   	
   xy<-matrix(csr(poly,1),1,2)
   delsq<-delta*delta
   while (dim(xy)[1]<n) {
      dsq<-0
      while (dsq<delsq) {
         xy.try<-c(csr(poly,1))
         dsq<-min((xy[,1]-xy.try[1])^2+(xy[,2]-xy.try[2])^2)
         }
      xy<-rbind(xy,xy.try)
      }
   if (k>0) {
   	  if(k > n/2) stop("k must be between 0 and n/2.")
      take<-matrix(sample(1:n,2*k,replace=FALSE),k,2)
      for (j in 1:k) {
         take1<-take[j,1]; take2<-take[j,2]
         xy1<-c(xy[take1,])
         angle<-2*pi*runif(1)
         radius<-rho*sqrt(runif(1))
         xy[take2,]<-xy1+radius*c(cos(angle),sin(angle))
         }
      }
    xy
} 

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom truncnorm rtruncnorm
##' @importFrom geoR varcov.spatial
binary.geo.Bayes <- function(formula,coords,data,
                               ID.coords,
                               control.prior,
                               control.mcmc,
                               kappa,messages) {
   kappa <- as.numeric(kappa)
   if(length(control.prior$log.prior.nugget)!=0 |
      length(control.mcmc$start.nugget) != 0) {
      stop("inclusion of the nugget effect is not available for binary data.")
   }
                                                                                               
   coords <- as.matrix(model.frame(coords,data,na.action=na.fail))

   out <- list()
   
   coords <- as.matrix(unique(coords))
   
   if(nrow(coords)!=nrow(unique(coords))) {
      warning("coincident locations may cause numerical issues in the computation: see ?jitterDupCoords")
   }

   mf <- model.frame(formula,data=data,na.action=na.fail)
   y <- as.numeric(model.response(mf))
   if(any(y > 1)) stop("the response variable must be binary")      
   n <- length(y)
   n.x <- nrow(coords)
   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))

   p <- ncol(D)
   
   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   
   U <- dist(coords)
   V.inv <- control.prior$V.inv
   m <- control.prior$m   	  
   nu <- 2*kappa
   n.S <- as.numeric(tapply(ID.coords,ID.coords,length))

   ind0 <- which(y==0)
   lim.y <- matrix(NA,ncol=2,nrow=n)
   lim.y[ind0,1] <- -Inf
   lim.y[ind0,2] <- 0
   lim.y[-ind0,1] <- 0
   lim.y[-ind0,2] <- Inf

   ind.beta <- 1:p
   ind.S <- (p+1):(n.x+p)
   cond.cov.inv <- matrix(NA,nrow=n.x+p,ncol=n.x+p)
   cond.cov.inv[ind.beta,ind.beta] <- V.inv+t(D)%*%D
   cond.cov.inv[ind.S,ind.beta] <- sapply(1:p,
                                   function(i) 
                                   tapply(D[,i],ID.coords,sum))
   cond.cov.inv[ind.beta,ind.S] <- t(cond.cov.inv[ind.S,ind.beta])
   beta.aux <- as.numeric(V.inv%*%m)

   full.cond.sim <- function(W,sigma2,R.inv) {
      cond.cov.inv[ind.S,ind.S] <- R.inv/sigma2 
      diag(cond.cov.inv[ind.S,ind.S]) <- 
      diag(cond.cov.inv[ind.S,ind.S])+n.S
      
      out.param.sim <- list()
      cond.cov <- solve(cond.cov.inv)
      b <- c(beta.aux+t(D)%*%W,tapply(W,ID.coords,sum))
      cond.mean <- as.numeric(cond.cov%*%b)
      out.sim <- list()
      sim.par <- cond.mean+t(chol(cond.cov))%*%rnorm(n.x+p)
      out.sim$beta <- sim.par[ind.beta]
      out.sim$S <- sim.par[ind.S]
      return(out.sim)
   }

   log.prior.theta1_2 <- function(theta1,theta2) {
      sigma2 <- exp(2*theta1)
      phi <- (exp(2*theta1-theta2))^(1/nu)
      log(phi)+log(sigma2)+control.prior$log.prior.sigma2(sigma2)+
      control.prior$log.prior.phi(phi)	
   }

   log.posterior <- function(theta1,theta2,S,param.post) {
	  sigma2 <- exp(2*theta1)
	  as.numeric(log.prior.theta1_2(theta1,theta2)+
	  -0.5*(n.x*log(sigma2)+param.post$ldetR+
	  t(S)%*%param.post$R.inv%*%(S)/sigma2))
   }

   acc.theta1 <- 0
   acc.theta2 <- 0

   n.sim <- control.mcmc$n.sim 
   burnin <- control.mcmc$burnin
   thin <- control.mcmc$thin
   h.theta1 <- control.mcmc$h.theta1
   h.theta2 <- control.mcmc$h.theta2
   
      
   c1.h.theta1 <- control.mcmc$c1.h.theta1
   c2.h.theta1 <- control.mcmc$c2.h.theta1
   c1.h.theta2 <- control.mcmc$c1.h.theta2
   c2.h.theta2 <- control.mcmc$c2.h.theta2
      
   
   out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
   colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")    

   out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=n.x)

   # Initialise all the parameters
   sigma2.curr <- control.mcmc$start.sigma2
   phi.curr <- control.mcmc$start.phi
	
   theta1.curr <- 0.5*log(sigma2.curr)
   theta2.curr <- log(sigma2.curr/(phi.curr^nu))
	
   beta.curr <- control.mcmc$start.beta
   mu.curr <- as.numeric(D%*%beta.curr)
   S.curr <- rep(0,n.x)
   W.curr <- rtruncnorm(n,a=lim.y[,1],b=lim.y[,2],
             mean=S.curr[ID.coords]+mu.curr)

   R.curr <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	         cov.pars=c(1,phi.curr),kappa=kappa,
	         nugget=0)$varcov
	  
  param.post.curr <- list()
  param.post.curr$R.inv <- solve(R.curr)
  param.post.curr$ldetR <- determinant(R.curr)$modulus
  lp.curr <- log.posterior(theta1.curr,theta2.curr,S.curr,param.post.curr)
     
       
  h1 <- rep(NA,n.sim) 
  h2 <- rep(NA,n.sim) 
      
  for(i in 1:n.sim) {
	
     # Update theta1 
     theta1.prop <- theta1.curr+h.theta1*rnorm(1)
     sigma2.prop <- exp(2*theta1.prop)
	 phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
     R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                          cov.pars=c(1,phi.prop),
	                          kappa=kappa,nugget=0)$varcov 
       
     param.post.prop <- list()	     
	 param.post.prop$R.inv <- solve(R.prop)
	 param.post.prop$ldetR <- determinant(R.prop)$modulus
	 lp.prop <- log.posterior(theta1.prop,theta2.curr,S.curr,param.post.prop) 
	 
	 if(log(runif(1)) < lp.prop-lp.curr) {
	    theta1.curr <- theta1.prop
	    param.post.curr <- param.post.prop
	    lp.curr <- lp.prop
	    sigma2.curr <- sigma2.prop
	    phi.curr <- phi.prop
	    acc.theta1 <- acc.theta1+1
	 } 
	 rm(theta1.prop,sigma2.prop,phi.prop)   	                                          
     h1[i] <- h.theta1 <- max(0,h.theta1 + 
	          (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                  
    	     
     # Update theta2 
     theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	 phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
     R.prop <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                          cov.pars=c(1,phi.prop),
	                          kappa=kappa,nugget=0)$varcov 
         
     param.post.prop <- list()	     
	 param.post.prop$R.inv <- solve(R.prop)
	 param.post.prop$ldetR <- determinant(R.prop)$modulus
	 lp.prop <- log.posterior(theta1.curr,theta2.prop,S.curr,
	                          param.post.prop) 	                                         
	     
	 if(log(runif(1)) < lp.prop-lp.curr) {
	    theta2.curr <- theta2.prop
	    param.post.curr <- param.post.prop
	    lp.curr <- lp.prop
	    phi.curr <- phi.prop
	    acc.theta2 <- acc.theta2+1
	 }    	                      
	 rm(theta2.prop,phi.prop)                     
     h2[i] <- h.theta2 <- max(0,h.theta2 + 
	          (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))  
	                                            
                    
     # Update beta, S and W
     sim.curr <- full.cond.sim(W.curr,sigma2.curr,param.post.curr$R.inv)   
     beta.curr <- sim.curr$beta
     mu.curr <- as.numeric(D%*%beta.curr)
     S.curr <- sim.curr$S
     W.curr <- rtruncnorm(n,lim.y[,1],lim.y[,2],
               mean=S.curr[ID.coords]+mu.curr)
     
     lp.curr <- log.posterior(theta1.curr,theta2.curr,S.curr,
	                          param.post.curr) 
	                           
     if(i > burnin & (i-burnin)%%thin==0) {
	    out$S[(i-burnin)/thin,] <- S.curr    
		out$estimate[(i-burnin)/thin,] <- c(beta.curr,sigma2.curr,
		                                    phi.curr)
	
   	 }
     if(messages) cat("Iteration",i,"out of",n.sim,"\r")
     flush.console()                            
        	
   }
   class(out) <- "Bayes.PrevMap"	
   out$y <- y
   out$D <- D
   out$coords <- coords
   out$kappa <- kappa
   out$ID.coords <- ID.coords  
   out$h1 <- h1
   out$h2 <- h2
   return(out)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom truncnorm rtruncnorm
##' @importFrom pdist pdist
binary.geo.Bayes.lr <- function(formula,coords,data,ID.coords,knots,
                       control.mcmc,control.prior,kappa,messages) {	
   knots <- as.matrix(knots)
   coords <- unique(as.matrix(model.frame(coords,data)))
   
      
   if(is.matrix(coords)==FALSE || dim(coords)[2] !=2) stop("coordinates must consist of two sets of numeric values.")                                                       	                                                       	
   mf <- model.frame(formula,data=data)
   y <- as.numeric(model.response(mf))
   n <- length(y)
   N <- nrow(knots)
   D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))	   
   p <- ncol(D)
   if(length(control.mcmc$start.beta)!=p) stop("wrong starting values for beta.")
   U.k <- as.matrix(pdist(coords,knots))
   V.inv <- control.prior$V.inv
   m <- control.prior$m   	  
   nu <- 2*kappa

   log.prior.theta1_2 <- function(theta1,theta2) {
      sigma2 <- exp(2*theta1)
      phi <- (exp(2*theta1-theta2))^(1/nu)
      log(phi/kappa)+control.prior$log.prior.sigma2(sigma2)+
      control.prior$log.prior.phi(phi)	
   }
             
   log.posterior <- function(theta1,theta2,mu,S,W,Z) {
	  sigma2 <- exp(2*theta1)
	  mu.Z <- as.numeric(mu+W[ID.coords])
	  diff.Z <- Z-mu.Z
	  as.numeric(log.prior.theta1_2(theta1,theta2)+
	  -0.5*(N*log(sigma2)+sum(S^2)/sigma2)+
     -0.5*sum(diff.Z^2))
   }
            
   c1.h.theta1 <- control.mcmc$c1.h.theta1
   c1.h.theta2 <- control.mcmc$c1.h.theta2
       
   c2.h.theta1 <- control.mcmc$c2.h.theta1
   c2.h.theta2 <- control.mcmc$c2.h.theta2
       
   h.theta1 <- control.mcmc$h.theta1
   h.theta2 <- control.mcmc$h.theta2
       
   acc.theta1 <- 0
   acc.theta2 <- 0
       
   n.sim <- control.mcmc$n.sim 
   burnin <- control.mcmc$burnin
   thin <- control.mcmc$thin
   
   out <- list()
   out$estimate <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=p+2)	
   colnames(out$estimate) <- c(colnames(D),"sigma^2","phi")
   	  	 
   out$S <- matrix(NA,nrow=(n.sim-burnin)/thin,ncol=N)
   out$const.sigma2 <- rep(NA,(n.sim-burnin)/thin)
       
   # Initialise all the parameters
   sigma2.curr <- control.mcmc$start.sigma2
   phi.curr <- control.mcmc$start.phi
   rho.curr <- phi.curr*2*sqrt(kappa)
   theta1.curr <- 0.5*log(sigma2.curr)
   theta2.curr <- log(sigma2.curr/(phi.curr^nu))

   beta.curr <- control.mcmc$start.beta
   mu.curr <- as.numeric(D%*%beta.curr)
   S.curr <- rep(0,N)

   ind0 <- which(y==0)
   lim.y <- matrix(NA,ncol=2,nrow=n)
   lim.y[ind0,1] <- -Inf
   lim.y[ind0,2] <- 0
   lim.y[-ind0,1] <- 0
   lim.y[-ind0,2] <- Inf

   ind.beta <- 1:p
   ind.S <- (p+1):(N+p)
   cond.cov.inv <- matrix(NA,nrow=N+p,ncol=N+p)
   cond.cov.inv[ind.beta,ind.beta] <- V.inv+t(D)%*%D
   beta.aux <- as.numeric(V.inv%*%m)

   n.S <- as.numeric(tapply(ID.coords,ID.coords,length))   
   D.aux <- sapply(1:p,function(i) tapply(D[,i],ID.coords,sum))

   full.cond.sim <- function(Z,sigma2,K) {
      cond.cov.inv[ind.S,ind.S] <- t(K*n.S)%*%K 
      diag(cond.cov.inv[ind.S,ind.S]) <- 
      diag(cond.cov.inv[ind.S,ind.S])+1/sigma2

      cond.cov.inv[ind.S,ind.beta] <- t(K)%*%D.aux
      cond.cov.inv[ind.beta,ind.S] <- t(cond.cov.inv[ind.S,ind.beta])

      out.param.sim <- list()
      cond.cov <- solve(cond.cov.inv)
      b <- c(beta.aux+t(D)%*%Z,t(K)%*%tapply(Z,ID.coords,sum))
      cond.mean <- as.numeric(cond.cov%*%b)
      out.sim <- list()
      sim.par <- cond.mean+t(chol(cond.cov))%*%rnorm(N+p)
      out.sim$beta <- sim.par[ind.beta]
      out.sim$S <- sim.par[ind.S]
      return(out.sim)
   }
   	
   # Compute the log-posterior density
   K.curr <- matern.kernel(U.k,rho.curr,kappa)	 
   W.curr <- as.numeric(K.curr%*%S.curr)
   Z.curr <- rtruncnorm(n,a=lim.y[,1],b=lim.y[,2],
             mean=mu.curr+W.curr[ID.coords],sd=1)
   lp.curr <- log.posterior(theta1.curr,theta2.curr,mu.curr,
              S.curr,W.curr,Z.curr)
      
       
   h1 <- rep(NA,n.sim)
   h2 <- rep(NA,n.sim)
   for(i in 1:n.sim) {

      # Update theta1 
      theta1.prop <- theta1.curr+h.theta1*rnorm(1)
      sigma2.prop <- exp(2*theta1.prop)
	  phi.prop <- (exp(2*theta1.prop-theta2.curr))^(1/nu)
	  rho.prop <- phi.prop*2*sqrt(kappa)
      K.prop <- matern.kernel(U.k,rho.curr,kappa)
         
      W.prop <- as.numeric(K.prop%*%S.curr)
	  lp.prop <- log.posterior(theta1.prop,theta2.curr,
	             mu.curr,S.curr,W.prop,Z.curr) 	                                         	     
	  if(log(runif(1)) < lp.prop-lp.curr) {
	     theta1.curr <- theta1.prop
	     lp.curr <- lp.prop
	     sigma2.curr <- sigma2.prop
	     phi.curr <- phi.prop
	     rho.curr <- rho.prop
	     K.curr <- K.prop
	     W.curr <- W.prop
	     acc.theta1 <- acc.theta1+1
	  } 
	  rm(theta1.prop,sigma2.prop,phi.prop,K.prop,W.prop,rho.prop)   	                                          
      h1[i] <- h.theta1 <- max(0,h.theta1 + 
	  (c1.h.theta1*i^(-c2.h.theta1))*(acc.theta1/i-0.45))                                    
    	     
      # Update theta2 
      theta2.prop <- theta2.curr+h.theta2*rnorm(1)
	  phi.prop <- (exp(2*theta1.curr-theta2.prop))^(1/nu)
	  rho.prop <- 2*phi.prop*sqrt(kappa)
      K.prop <- matern.kernel(U.k,rho.prop,kappa) 
         
      W.prop <- as.numeric(K.prop%*%S.curr)
	  lp.prop <- log.posterior(theta1.curr,theta2.prop,
	             mu.curr,S.curr,W.prop,Z.curr) 	                                         
	     
	  if(log(runif(1)) < lp.prop-lp.curr) {
	     theta2.curr <- theta2.prop
	     lp.curr <- lp.prop
	     phi.curr <- phi.prop
	     rho.curr <- rho.prop
	     K.curr <- K.prop
	     W.curr <- W.prop
	     acc.theta2 <- acc.theta2+1
	  }    	                      
	  rm(theta2.prop,phi.prop,rho.prop,K.prop,W.prop)                     
      h2[i] <- h.theta2 <- max(0,h.theta2 + 
	  (c1.h.theta2*i^(-c2.h.theta2))*(acc.theta2/i-0.45))                                                          
	                          
      # Update beta, S and Z
      sim.curr <- full.cond.sim(Z.curr,sigma2.curr,K.curr)
      beta.curr <- sim.curr$beta
      mu.curr <- as.numeric(D%*%beta.curr)
      S.curr <- sim.curr$S
      W.curr <- K.curr%*%S.curr
      Z.curr <- rtruncnorm(n,a=lim.y[,1],b=lim.y[,2],
                mean=mu.curr+W.curr[ID.coords],sd=1)
      
      lp.curr <- log.posterior(theta1.curr,theta2.curr,
	             mu.curr,S.curr,W.curr,Z.curr)   
	                           
      if(i > burnin & (i-burnin)%%thin==0) {
	     out$S[(i-burnin)/thin,] <- S.curr
		 out$estimate[(i-burnin)/thin,] <-
		 c(beta.curr,sigma2.curr,phi.curr)	
		 out$const.sigma2[(i-burnin)/thin] <- 
		 mean(apply(K.curr,1,function(r) sqrt(sum(r^2))))
   	  }
      if(messages) cat("Iteration",i,"out of",n.sim,"\r")
      flush.console()                            
   }

   class(out) <- "Bayes.PrevMap"	
   out$y <- y
   out$D <- D
   out$coords <- coords       
   out$kappa <- kappa  
   out$knots <- knots   
   out$h1 <- h1
   out$h2 <- h2
   return(out)
}


##' @title Bayesian estimation for the two-levels binary probit model
##' @description This function performs Bayesian estimation for a geostatistical binary probit model. It also allows to specify a two-levels model so as to include individual-level and household-level (or any other unit comprising a group of individuals, e.g. village, school, compound, etc...) variables. 
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param ID.coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID.coords}}. These must be provided in order to specify spatial random effects at household-level. \bold{Warning}: the household coordinates must all be distinct otherwise see \code{\link{jitterDupCoords}}. Default is \code{NULL}.
##' @param control.prior output from \code{\link{control.prior}}.
##' @param control.mcmc output from \code{\link{control.mcmc.Bayes}}.
##' @param kappa value for the shape parameter of the Matern covariance function.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation is required. Default is \code{low.rank=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @details
##' This function performs Bayesian estimation for the parameters of the geostatistical binary probit model. Let \eqn{i} and \eqn{j} denote the indices of the \eqn{i}-th household and \eqn{j}-th individual within that household. The response variable \eqn{Y_{ij}} is a binary indicator taking value 1 if the individual has been tested positive for the disease of interest and 0 otherwise. Conditionally on a zero-mean stationary Gaussian process \eqn{S(x_{i})}, \eqn{Y_{ij}} are mutually independent Bernoulli variables with probit link function \eqn{\Phi^{-1}(\cdot)}, i.e.
##' \deqn{\Phi^{-1}(p_{ij}) = d_{ij}'\beta + S(x_{i}),}
##' where \eqn{d_{ij}} is a vector of covariates, both at individual- and household-level, with associated regression coefficients \eqn{\beta}. The Gaussian process \eqn{S(x)} has isotropic Matern covariance function (see \code{\link{matern}}) with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. 
##' 
##' \bold{Priors definition.} Priors can be defined through the function \code{\link{control.prior}}. The hierarchical structure of the priors is the following. Let \eqn{\theta} be the vector of the covariance parameters \code{c(sigma2,phi)}; each component of \eqn{\theta} has independent priors that can be freely defined by the user. However, in  \code{\link{control.prior}} uniform and log-normal priors are also available as default priors for each of the covariance parameters. The vector of regression coefficients \code{beta} has a multivariate Gaussian prior with mean \code{beta.mean} and covariance matrix \code{beta.covar}. 
##'
##' \bold{Updating regression coefficents and random effects using auxiliary variables.} To update \eqn{\beta} and \eqn{S(x_{i})}, we use an auxiliary variable technique based on Rue and Held (2005). Let \eqn{V_{ij}} denote a set of random variables that conditionally on \eqn{\beta} and \eqn{S(x_{i})}, are mutually independent Gaussian with mean \eqn{d_{ij}'\beta + S(x_{i})} and unit variance. Then, \eqn{Y_{ij}=1} if \eqn{V_{ij} > 0} and \eqn{Y_{ij}=0} otherwise. Using this representation of the model, we use a Gibbs sampler to simulate from the full conditionals of \eqn{\beta}, \eqn{S(x_{i})} and \eqn{V_{ij}}. See Section 4.3 of Rue and Held (2005) for more details. 
##'
##' \bold{Updating the covariance parameters with a Metropolis-Hastings algorithm.} In the MCMC algorithm implemented in \code{binary.probit.Bayes}, the transformed parameters \deqn{(\theta_{1}, \theta_{2})=(\log(\sigma^2)/2,\log(\sigma^2/\phi^{2 \kappa}))} are independently updated using a Metropolis Hastings algorithm. At the \eqn{i}-th iteration, a new value is proposed for each parameter from a univariate Gaussian distrubion with variance \eqn{h_{i}^2}. This is tuned using the following adaptive scheme \deqn{h_{i} = h_{i-1}+c_{1}i^{-c_{2}}(\alpha_{i}-0.45),} where \eqn{\alpha_{i}} is the acceptance rate at the \eqn{i}-th iteration, 0.45 is the optimal acceptance rate for a univariate Gaussian distribution, whilst \eqn{c_{1} > 0} and \eqn{0 < c_{2} < 1} are pre-defined constants. The starting values \eqn{h_{1}} for each of the parameters \eqn{\theta_{1}} and \eqn{\theta_{2}} can be set using the function \code{\link{control.mcmc.Bayes}} through the arguments \code{h.theta1}, \code{h.theta2} and \code{h.theta3}. To define values for \eqn{c_{1}} and \eqn{c_{2}}, see the documentation of \code{\link{control.mcmc.Bayes}}.
##' 
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), \code{sigma2} may take very different values from the actual variance of the Gaussian process to approximate. The function \code{\link{adjust.sigma2}} can then be used to (approximately) explore the range for \code{sigma2}. For example if the variance of the Gaussian process is \code{0.5}, then an approximate value for \code{sigma2} is \code{0.5/const.sigma2}, where \code{const.sigma2} is the value obtained with \code{\link{adjust.sigma2}}.
##' @return An object of class "Bayes.PrevMap".
##' The function \code{\link{summary.Bayes.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: matrix of the posterior samples of the model parameters.
##' @return \code{S}: matrix of the posterior samples for each component of the random effect. 
##' @return \code{const.sigma2}: vector of the values of the multiplicative factor used to adjust the values of \code{sigma2} in the low-rank approximation.
##' @return \code{y}: binary observations.
##' @return \code{D}: matrix of covariarates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{kappa}: shape parameter of the Matern function.
##' @return \code{ID.coords}: set of ID values defined through the argument \code{ID.coords}.
##' @return \code{knots}: matrix of spatial knots used in the low-rank approximation.
##' @return \code{h1}: vector of values taken by the tuning parameter \code{h.theta1} at each iteration.
##' @return \code{h2}: vector of values taken by the tuning parameter \code{h.theta2} at each iteration.
##' @return \code{call}: the matched call.
##' @seealso  \code{\link{control.mcmc.Bayes}},  \code{\link{control.prior}},\code{\link{summary.Bayes.PrevMap}}, \code{\link{matern}}, \code{\link{matern.kernel}}, \code{\link{create.ID.coords}}.
##' @references Rue, H., Held, L. (2005). \emph{Gaussian Markov Random Fields: Theory and Applications.} Chapman & Hall, London.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
binary.probit.Bayes <- function(formula,coords,data,ID.coords,
                           control.prior,control.mcmc,kappa,low.rank=FALSE,
                           knots=NULL,messages=TRUE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.") 
    if(class(control.mcmc)!="mcmc.Bayes.PrevMap") stop("control.mcmc must be of class 'mcmc.Bayes.PrevMap'")

    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(kappa < 0) stop("kappa must be positive.")
	if(!low.rank) {
		res <- binary.geo.Bayes(formula=formula,coords=coords,
		       data=data,ID.coords=ID.coords,
		       control.prior=control.prior,
		       control.mcmc=control.mcmc,
		       kappa=kappa,messages=messages)
	} else {
		res <- binary.geo.Bayes.lr(formula=formula,
		       coords=coords,data=data,ID.coords=ID.coords,
		       knots=knots,control.mcmc=control.mcmc,
		       control.prior=control.prior,kappa=kappa,
		       messages=messages)
	}
	res$call <- match.call()
	return(res)
}

##' @title Monte Carlo Maximum Likelihood estimation for the Poisson model
##' @description This function performs Monte Carlo maximum likelihood (MCML) estimation for the geostatistical Poisson model with log link function.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param units.m an object of class \code{\link{formula}} indicating the multiplicative offset for the mean of the Poisson model; if not specified this is then internally set as 1.
##' @param coords an object of class \code{\link{formula}} indicating the geographic coordinates.
##' @param data a data frame containing the variables in the model. 
##' @param par0 parameters of the importance sampling distribution: these should be given in the following order \code{c(beta,sigma2,phi,tau2)}, where \code{beta} are the regression coefficients, \code{sigma2} is the variance of the Gaussian process, \code{phi} is the scale parameter of the spatial correlation and \code{tau2} is the variance of the nugget effect (if included in the model).
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param kappa fixed value for the shape parameter of the Matern covariance function.
##' @param fixed.rel.nugget fixed value for the relative variance of the nugget effect; \code{fixed.rel.nugget=NULL} if this should be included in the estimation. Default is \code{fixed.rel.nugget=NULL}.
##' @param start.cov.pars a vector of length two with elements corresponding to the starting values of \code{phi} and the relative variance of the nugget effect \code{nu2}, respectively, that are used in the optimization algorithm. If \code{nu2} is fixed through \code{fixed.rel.nugget}, then \code{start.cov.pars} represents the starting value for \code{phi} only.
##' @param method method of optimization. If \code{method="BFGS"} then the \code{\link{maxBFGS}} function is used; otherwise \code{method="nlminb"} to use the \code{\link{nlminb}} function. Default is \code{method="BFGS"}.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation of the Gaussian spatial process is used when fitting the model. Default is \code{low.rank=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots that are used in the low-rank approximation. Default is \code{knots=NULL}. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the samples of the random effect is displayed after completion of conditional simulation. Default is \code{plot.correlogram=TRUE}.
##' @details
##' This function performs parameter estimation for a geostatistical Poisson model with log link function. Conditionally on a zero-mean stationary Gaussian process \eqn{S(x)} and mutually independent zero-mean Gaussian variables \eqn{Z} with variance \code{tau2}, the observations \code{y} are generated from a Poisson distribution with mean \eqn{m\lambda}, where \eqn{m} is an offset defined through the argument \code{units.m}. A canonical log link is used, thus the linear predictor assumes the form
##' \deqn{\log(\lambda) = d'\beta + S(x) + Z,}
##' where \eqn{d} is a vector of covariates with associated regression coefficients \eqn{\beta}. The Gaussian process \eqn{S(x)} has isotropic Matern covariance function (see \code{\link{matern}}) with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}. 
##' In the \code{poisson.log.MCML} function, the shape parameter is treated as fixed. The relative variance of the nugget effect, \code{nu2=tau2/sigma2}, can also be fixed through the argument \code{fixed.rel.nugget}; if \code{fixed.rel.nugget=NULL}, then the relative variance of the nugget effect is also included in the estimation.
##' 
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)+Z} given the data \code{y}, in order to approximate the high-dimensiional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic epression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{maxBFGS}} (\code{method="BFGS"}), from the \pkg{maxLik} package, and \code{\link{nlminb}} (\code{method="nlminb"}).
##' 
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), the parameter \code{sigma2} is then multiplied by a factor \code{constant.sigma2} so as to obtain a value that is closer to the actual variance of \eqn{S(x)}. 
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: estimates of the model parameters; use the function \code{\link{coef.PrevMap}} to obtain estimates of covariance parameters on the original scale.
##' @return \code{covariance}: covariance matrix of the MCML estimates.
##' @return \code{log.lik}: maximum value of the log-likelihood.
##' @return \code{y}: observations.
##' @return \code{units.m}: offset.
##' @return \code{D}: matrix of covariates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{method}: method of optimization used.
##' @return \code{kappa}: fixed value of the shape parameter of the Matern function.
##' @return \code{knots}: matrix of the spatial knots used in the low-rank approximation.
##' @return \code{const.sigma2}: adjustment factor for \code{sigma2} in the low-rank approximation.
##' @return \code{h}: vector of the values of the tuning parameter at each iteration of the Langevin-Hastings MCMC algorithm; see \code{\link{Laplace.sampling}}, or \code{\link{Laplace.sampling.lr}} if a low-rank approximation is used.
##' @return \code{samples}: matrix of the random effects samples from the importance sampling distribution used to approximate the likelihood function.
##' @return \code{fixed.rel.nugget}: fixed value for the relative variance of the nugget effect.
##' @return \code{call}: the matched call.
##' @seealso \code{\link{Laplace.sampling}}, \code{\link{Laplace.sampling.lr}}, \code{\link{summary.PrevMap}}, \code{\link{coef.PrevMap}}, \code{\link{matern}}, \code{\link{matern.kernel}},  \code{\link{control.mcmc.MCML}}.
##' @references Christensen, O. F. (2004). \emph{Monte carlo maximum likelihood in model-based geostatistics.} Journal of Computational and Graphical Statistics 13, 702-718.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
poisson.log.MCML <- function(formula,units.m=NULL,coords,data,
                             par0,control.mcmc,kappa,
                             fixed.rel.nugget=NULL,
                             start.cov.pars,
                             method="BFGS",low.rank=FALSE,
                             knots=NULL,
                             messages=TRUE,
                             plot.correlogram=TRUE) {
    if(low.rank & length(dim(knots))==0) stop("if low.rank=TRUE, then knots must be provided.") 
    if(class(control.mcmc)!="mcmc.MCML.PrevMap") stop("control.mcmc must be of class 'mcmc.MCML.PrevMap'")
    if(class(formula)!="formula") stop("formula must be a 'formula' object indicating the variables of the model to be fitted.")
    if(class(coords)!="formula") stop("coords must be a 'formula' object indicating the spatial coordinates.")
    if(length(units.m)>0 && class(units.m)!="formula") stop("units.m must be a 'formula' object indicating the offset for the mean of the Poisson model.")
    if(kappa < 0) stop("kappa must be positive.")
    if(method != "BFGS" & method != "nlminb") stop("'method' must be either 'BFGS' or 'nlminb'.")
    if(!low.rank) {
		res <- geo.MCML(formula=formula,units.m=units.m,coords=coords,
		           data=data,ID.coords=NULL,par0=par0,control.mcmc=control.mcmc,
		           kappa=kappa,fixed.rel.nugget=fixed.rel.nugget,
		           start.cov.pars=start.cov.pars,method=method,messages=messages,
		           plot.correlogram=plot.correlogram,poisson.llik=TRUE)
    } else {
		res <- geo.MCML.lr(formula=formula,units.m=units.m,coords=coords,
		           data=data,knots=knots,par0=par0,control.mcmc=control.mcmc,
		           kappa=kappa,start.cov.pars=start.cov.pars,method=method,
		           messages=messages,plot.correlogram=plot.correlogram,poisson.llik=TRUE)		
    } 
    res$call <- match.call()
    return(res)
}


##' @title Spatial predictions for the Poisson model with log link function, using plug-in of MCML estimates
##' @description This function performs spatial prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a geostatistical Poisson model with log link function.
##' @param object an object of class "PrevMap" obtained as result of a call to \code{\link{poisson.log.MCML}}.
##' @param grid.pred a matrix of prediction locations.
##' @param predictors a data frame of the values of the explanatory variables at each of the locations in \code{grid.pred}; each column correspond to a variable and each row to a location. \bold{Warning:} the names of the columns in the data frame must match those in the data used to fit the model. Default is \code{predictors=NULL} for models with only an intercept.
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param type a character indicating the type of spatial predictions: \code{type="marginal"} for marginal predictions or \code{type="joint"} for joint predictions. Default is \code{type="marginal"}. In the case of a low-rank approximation only joint predictions are available.
##' @param scale.predictions a character vector of maximum length 2, indicating the required scale on which spatial prediction is carried out: "log" and "exponential". Default is \code{scale.predictions=c("log","exponential")}.
##' @param quantiles a vector of quantiles used to summarise the spatial predictions.
##' @param standard.errors logical; if \code{standard.errors=TRUE}, then standard errors for each \code{scale.predictions} are returned. Default is \code{standard.errors=FALSE}.
##' @param thresholds a vector of exceedance thresholds; default is \code{thresholds=NULL}.
##' @param scale.thresholds a character value indicating the scale on which exceedance thresholds are provided; \code{"log"} or \code{"exponential"}. Default is \code{scale.thresholds=NULL}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the conditional simulations is displayed. 
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @return A "pred.PrevMap" object list with the following components: \code{log}; \code{exponential}; \code{exceedance.prob}, corresponding to a matrix of the exceedance probabilities where each column corresponds to a specified value in \code{thresholds}; \code{samples}, corresponding to a matrix of the predictive samples at each prediction locations for the linear predictor of the Poisson model (if \code{scale.predictions="log"} this component is \code{NULL}); \code{grid.pred} prediction locations. 
##' Each of the three components \code{log} and  \code{exponential} is also a list with the following components:
##' @return \code{predictions}: a vector of the predictive mean for the associated quantity (log or exponential).
##' @return \code{standard.errors}: a vector of prediction standard errors (if \code{standard.errors=TRUE}).
##' @return \code{quantiles}: a matrix of quantiles of the resulting predictions with each column corresponding to a quantile specified through the argument \code{quantiles}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk} 
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom geoR matern 
##' @export 

spatial.pred.poisson.MCML <- function(object,grid.pred,predictors=NULL,control.mcmc,
                                     type="marginal",
                                     scale.predictions=c("log","exponential"),
                                     quantiles=c(0.025,0.975),
                                     standard.errors=FALSE,                                                                                                                                                      
                                     thresholds=NULL,
                                     scale.thresholds=NULL,
                                     plot.correlogram=FALSE,
                                     messages=TRUE) {
    if(nrow(grid.pred) < 2) stop("prediction locations must be at least two.")       
    if(length(predictors)>0 && class(predictors)!="data.frame") stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.")
    if(length(predictors)>0 && any(is.na(predictors))) stop("missing values found in 'predictors'.")
    p <- object$p <- ncol(object$D)	
    kappa <- object$kappa	
    n.pred <- nrow(grid.pred)
    coords <- object$coords
    if(type=="marginal" & length(object$knots) > 0) warning("only joint predictions are avilable for the low-rank approximation.")
    if(any(type==c("marginal","joint"))==FALSE) stop("type of predictions should be marginal or joint")
    for(i in 1:length(scale.predictions)) {
       if(any(c("log","exponential")==scale.predictions[i])==FALSE) stop("invalid scale.predictions.")
    }
	
    if(length(thresholds)>0) {
       if(any(scale.predictions==scale.thresholds)==FALSE) {
	    stop("scale thresholds must be equal to a scale prediction")
       }
    }
	
    if(length(thresholds)==0 & length(scale.thresholds)>0 |
       length(thresholds)>0 & length(scale.thresholds)==0) stop("to estimate exceedance probabilities both thresholds and scale.thresholds.")

    if(object$p==1) {
       predictors <- matrix(1,nrow=n.pred)	
    } else {
	 if(length(dim(predictors))==0) stop("covariates at prediction locations should be provided.")	
	 predictors <- as.matrix(model.matrix(delete.response(terms(formula(object$call))),data=predictors))
	 if(nrow(predictors)!=nrow(grid.pred)) stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")
	 if(ncol(predictors)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
    }
	
    if(length(dim(object$knots)) > 0) {	
	 beta <- object$estimate[1:p]
	 sigma2 <- exp(object$estimate["log(sigma^2)"])/object$const.sigma2
	 rho <- exp(object$estimate["log(phi)"])*2*sqrt(object$kappa)
	 knots <- object$knots
	 U.k <- as.matrix(pdist(coords,knots))
	 K <- matern.kernel(U.k,rho,kappa)
	 mu.pred <- as.numeric(predictors%*%beta)
	 object$mu <- object$D%*%beta	
	 Z.sim.res <- Laplace.sampling.lr(object$mu,sigma2,K,
	 object$y,object$units.m,control.mcmc,
	 plot.correlogram=plot.correlogram,messages=messages,poisson.llik=TRUE)
	 Z.sim <- Z.sim.res$samples
	 U.k.pred <- as.matrix(pdist(grid.pred,knots))
	 K.pred <- matern.kernel(U.k.pred,rho,kappa) 
	 out <- list()     
       eta.sim <- sapply(1:(dim(Z.sim)[1]), function(i) mu.pred+K.pred%*%Z.sim[i,])           
    
       if(any(scale.predictions=="log")) {
    	    if(messages) cat("Spatial predictions: log \n")    	   
    	    out$log$predictions <- apply(eta.sim,1,mean)
          if(standard.errors) {
             out$log$standard.errors <- apply(eta.sim,1,sd)
          }
          if(length(quantiles) > 0) {
             out$log$quantiles <- t(apply(eta.sim,1,function(r) quantile(r,quantiles)))         
          }
       
          if(length(thresholds) > 0 && scale.thresholds=="log") {
             out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
             colnames(out$exceedance.prob) <- paste(thresholds,sep="")
             for(j in 1:length(thresholds)) {
               out$exceedance.prob[,j] <- apply(eta.sim,1,function(r) mean(r > thresholds[j]))	
             } 
          }		   
       }
       
       if(any(scale.predictions=="exponential")) {
          if(messages) cat("Spatial predictions: exponential \n") 
          exponential.sim <- exp(eta.sim)    	   
          out$exponential$predictions <- apply(exponential.sim,1,mean)
          if(standard.errors) {
             out$exponential$standard.errors <- apply(exponential.sim,1,sd)
          }
       
          if(length(quantiles) > 0) {
             out$exponential$quantiles <- t(apply(exponential.sim,1,function(r) quantile(r,quantiles)))         
          }
       	
          if(length(thresholds) > 0 && scale.thresholds=="exponential") {
             out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	 colnames(out$exceedance.prob) <- paste(thresholds,sep="")
             for(j in 1:length(thresholds)) {
                out$exceedance.prob[,j] <- apply(exponential.sim,1,function(r) mean(r > thresholds[j]))	
             } 
          }	   
       }  
    } else {	
       beta <- object$estimate[1:p]
	 sigma2 <- exp(object$estimate[p+1])
	 phi <- exp(object$estimate[p+2])
	 if(length(object$fixed.rel.nugget)==0){		
	    tau2 <- sigma2*exp(object$estimate[p+3]) 
	 } else {
	    tau2 <- object$fixed.rel.nugget*sigma2
	 }
	

	 U <- dist(coords)
	 U.pred.coords <- as.matrix(pdist(grid.pred,coords))	
	 Sigma <- varcov.spatial(dists.lowertri=U,cov.model="matern",
	                cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)$varcov 
	 Sigma.inv <- solve(Sigma)   	      
       C <- sigma2*matern(U.pred.coords,phi,kappa)	                
       A <- C%*%Sigma.inv
       out <- list()  
	 mu.pred <- as.numeric(predictors%*%beta)   
	 object$mu <- object$D%*%beta      	             	
    
	 S.sim.res <- Laplace.sampling(object$mu,Sigma,
	              object$y,object$units.m,control.mcmc,
	              plot.correlogram=plot.correlogram,messages=messages,poisson.llik=TRUE) 	   	
       S.sim <- S.sim.res$samples     
       mu.cond <- sapply(1:(dim(S.sim)[1]),function(i) mu.pred+A%*%(S.sim[i,]-object$mu))   	
    
	              
       if(type=="marginal") {
    	    if(messages) cat("Type of predictions:",type,"\n")
    	    sd.cond <- sqrt(sigma2-diag(A%*%t(C)))
       } else if (type=="joint") {
    	    if(messages) cat("Type of predictions: ",type," (this step might be demanding) \n")
    	    Sigma.pred <-  varcov.spatial(coords=grid.pred,cov.model="matern",
	                  cov.pars=c(sigma2,phi),kappa=kappa)$varcov 
	    Sigma.cond <- Sigma.pred - A%*%t(C)
	    sd.cond <- sqrt(diag(Sigma.cond))	   
       }  
    
       if((length(quantiles) > 0) | 
          (any(scale.predictions=="exponential")) |
          (length(scale.thresholds)>0)) {
          if(type=="marginal") {
             eta.sim <- sapply(1:(dim(S.sim)[1]), function(i) rnorm(n.pred,mu.cond[,i],sd.cond))
          } else if(type=="joint") {
             Sigma.cond.sroot <- t(chol(Sigma.cond))
             eta.sim <- sapply(1:(dim(S.sim)[1]), function(i) mu.cond[,i]+
                                                                   Sigma.cond.sroot%*%rnorm(n.pred))                                                                   
          }      	
       }
    
       if(any(scale.predictions=="log")) {
    	    if(messages) cat("Spatial predictions: log \n")    	   
    	       out$log$predictions <- apply(mu.cond,1,mean)
          if(standard.errors) {
             out$log$standard.errors <- sqrt(sd.cond^2+diag(A%*%cov(S.sim)%*%t(A)))
          }
          
          if(length(quantiles) > 0) {
             out$log$quantiles <- t(apply(eta.sim,1,function(r) quantile(r,quantiles)))         
          }
       
          if(length(thresholds) > 0 && scale.thresholds=="log") {
       	 out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	 colnames(out$exceedance.prob) <- paste(thresholds,sep="")
             for(j in 1:length(thresholds)) {
                out$exceedance.prob[,j] <- apply(eta.sim,1,function(r) mean(r > thresholds[j]))	
             } 
          }		   
       
       }
       
       if(any(scale.predictions=="exponential")) {
          if(messages) cat("Spatial predictions: exponential \n") 
          exponential.sim <- exp(eta.sim)    	   
          out$exponential$predictions <- apply(exp(mu.cond+0.5*sd.cond^2),1,mean)
          if(standard.errors) {
             out$exponential$standard.errors <- apply(exponential.sim,1,sd)
          }
       
          if(length(quantiles) > 0) {
             out$exponential$quantiles <- t(apply(exponential.sim,1,function(r) quantile(r,quantiles)))         
          }
       	
          if(length(thresholds) > 0 && scale.thresholds=="exponential") {
       	 out$exceedance.prob <- matrix(NA,nrow=n.pred,ncol=length(thresholds))
       	 colnames(out$exceedance.prob) <- paste(thresholds,sep="")
             for(j in 1:length(thresholds)) {
                out$exceedance.prob[,j] <- apply(exponential.sim,1,function(r) mean(r > thresholds[j]))	
             } 
          }	   
       }
    }
    
    if(any(scale.predictions=="exponential")) {
       out$samples <- eta.sim	
    }
    out$grid <- grid.pred
    class(out) <- "pred.PrevMap"
    out
}
