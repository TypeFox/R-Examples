pchisqsum2 <- function(Q, lambda, delta = rep(0,length(lambda)), method=c("saddlepoint","integration","liu"),acc=1e-7){
  method <- match.arg(method)
  
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
	
  if(method =="saddlepoint"){
  		p = saddle(Q,lambda,delta)
		if(is.na(p)){
			method <- "integration"
		} else { 
			return(list(p =p, errflag=0))
		}
	}
	
  if(method == "integration"){
		tmp <- davies(q=Q, lambda=lambda, delta=delta, acc=acc)
		if(tmp$ifault > 0){
		  lambda <- zapsmall(lambda, digits=2)
		  delta <- delta[lambda > 0]
		  lambda <- lambda[lambda > 0]
      	  tmp <- farebrother(q=Q,lambda=lambda,delta=delta)
		  tmp$Qq <- tmp$res
		}
    return(list(p = tmp$Qq, errflag = 0))
	}
	
  if(method == "liu"){
		tmp <- liu(Q, lambda=lambda, delta=delta )
		return(list(p = tmp, errflag = 0))
	}
}


saddle <- function (x, lambda,delta=rep(0,length(lambda))){
   	if(x==0) return(1)
    if(length(lambda) ==1){
      return(stats::pchisq(x/lambda,df=1,ncp=delta,lower.tail=FALSE))
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta){
    	-sum(log(1 - 2 * zeta * lambda))/2 + sum((delta*lambda*zeta)/(1-2*zeta*lambda))
    	}
    kprime0 <- function(zeta){
    	sapply(zeta, function(zz){ 
    		sum(lambda/(1 - 2 * zz * lambda)) + sum((delta*lambda)/(1-2*zeta*lambda)+2*(delta*zz*lambda^2)/(1-2*zz*lambda)^2)
    		})
    	}
    kpprime0 <- function(zeta){ 
    	2 * sum(lambda^2/(1 - 2 * zeta * lambda)^2) +sum((4*delta*lambda^2)/(1-2*zeta*lambda)^2+ 8*delta*zeta*lambda^3/(1-2*zeta*lambda)^3)
    	}
    n <- length(lambda)
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum(lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- stats::uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, 
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04) return(NA)
    else return(stats::pnorm(w + log(v/w)/w, lower.tail = FALSE))
}