##################################################################
######  "Random number generators and 
######          related density functions for yuima packages"
##################################################################


## Bilateral gamma

rbgamma <- function(x,delta.plus,gamma.plus,delta.minus,gamma.minus){
 if( delta.plus <= 0 )
   stop("delta.plus must be positive.")
 if( gamma.plus <= 0 )
   stop("gamma.plus must be positive.")
 if( delta.minus <= 0 )
   stop("delta.minus must be positive.")
 if( gamma.minus <= 0 )
   stop("gamma.minus must be positive.")
 X <- rgamma(x,delta.plus,gamma.plus) - rgamma(x,delta.minus,gamma.minus)
 return(X)
}

# "dbgamma" by YU.
dbgamma<-function(x,delta.plus,gamma.plus,delta.minus,gamma.minus){
  ## Error check
  if(length(delta.plus)!=1||length(gamma.plus)!=1||length(delta.minus)!=1||length(gamma.minus)!=1)
  stop("All of the parameters are numeric.")
  if(delta.plus<=0||gamma.plus<=0||delta.minus<=0||gamma.minus<=0)
  stop("All of the parameters are positive.")
  
  ## On the positive line
  funcp<-function(x,y){y^{delta.minus-1}*(x+y/(gamma.plus+gamma.minus))^{delta.plus-1}*exp(-y)}
  intp<-function(x){integrate(funcp,lower=0,upper=Inf,x=x)$value}
  intvecp<-function(x)sapply(x,intp)
  densp<-gamma.plus^delta.plus*gamma.minus^delta.minus/((gamma.plus+gamma.minus)^delta.minus*gamma(delta.plus)*gamma(delta.minus))*exp(-gamma.plus*x)*intvecp(x)
  
  ## On the negative line
  funcm<-function(x,y){y^{delta.plus-1}*(-x+y/(gamma.plus+gamma.minus))^{delta.minus-1}*exp(-y)}
  intm<-function(x){integrate(funcm,lower=0,upper=Inf,x=x)$value}
  intvecm<-function(x)sapply(x,intm)
  densm<-gamma.plus^delta.plus*gamma.minus^delta.minus/((gamma.plus+gamma.minus)^delta.plus*gamma(delta.minus)*gamma(delta.plus))*exp(gamma.minus*x)*intvecm(x)

  dens<-ifelse(0<=x,densp,densm)
  dens
  
}


## (Multivariate) Normal gamma

rngamma <- function(x,lambda,alpha,beta,mu,Lambda){
  ## Error check
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
   Lambda <- NA
  
  if(is.na(Lambda)){
    ## univariate case
    if(length(mu)!=1 || length(beta)!=1){
      stop("Error: wrong input dimension.")
    }
    tmp <- as.numeric(alpha^2 - beta^2)
    if( lambda <= 0 ){
      stop("lambda must be positive.")
    }
    if( alpha <= 0 ){
      stop("alpha must be positive.")
    }
    if( tmp <= 0 ){
      stop("alpha^2 - beta^2 must be positive value.")
    }
  
    tau <- rgamma(x,lambda,tmp/2)
    eta <- rnorm(x)
    ##  z <- mu + beta * tau * Lambda + sqrt(tau * Lambda) * eta
    z <- mu + beta * tau + sqrt(tau) * eta
    X <- z
    return(X)
    
  }else{ ## multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    
    if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }

    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
  
    if( lambda <= 0 )
      stop("lambda must be positive.")
    if( alpha <= 0 )
      stop("alpha must be positive.")
    if( tmp <=0)
      stop("alpha^2 - t(beta) %*% Lambda %*% must be positive.")
    
    tau <- rgamma(x,lambda,tmp/2)
    eta <- rnorm(x*length(beta))
    sqrt.L <- svd(Lambda)
    sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)
    
    z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
    X <- z
    return(X)
  }
}


dngamma <- function(x,lambda,alpha,beta,mu,Lambda){
  ## Error check
  if(length(lambda)!=1||length(alpha)!=1)
    stop("alpha and lambda must be positive reals.")
  if( lambda <= 0 )
    stop("lambda must be positive.")
  if( alpha <= 0 )
    stop("alpha must be positive.")
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
    Lambda <- NA
  if(is.na(Lambda)){
    ## univariate case
    if(length(mu)!=1 || length(beta)!=1){
      stop("Error: wrong input dimension.")
    }
    if( alpha^2 - beta^2 <= 0 )
      stop("alpha^2 - beta^2 must be positive.")
    
    dens <- exp(beta*(x-mu))*((alpha^2 - beta^2)^(lambda))*besselK(alpha*abs(x-mu),lambda-1/2)*abs(x-mu)^(lambda-1/2)/(gamma(lambda)*sqrt(pi)*((2*alpha)^(lambda-1/2)))
    dens}
  else{ ## multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    
    if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }

    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
    if( tmp <=0)
      stop("alpha^2 - t(beta) %*% Lambda %*% must be positive.")
    Lambdainv<-solve(Lambda)
    dens<- exp(t(beta)%*%(x-mu))*(alpha^2-t(beta)%*%Lambda%*%beta)^(lambda)*besselK(alpha*sqrt(t(x-mu)%*%Lambdainv%*%(x-mu)),lambda-nrow(Lambda)/2)*sqrt(t(x-mu)%*%Lambdainv%*%(x-mu))^{lambda-nrow(Lambda)/2}/(gamma(lambda)*pi^{nrow(Lambda)/2}*2^{nrow(Lambda)/2+lambda-1}*alpha^{lambda-nrow(Lambda)/2})
    dens
  }
}


## Inverse Gaussian

rIG <- function(x,delta,gamma){
  if( delta <= 0 )
    stop("delta must be positive.")
  if( gamma <= 0 )
    stop("gamma must be positive.")
  V <- rchisq(x,df=1);
  
  z1 <- ( delta/gamma + V/(2*gamma^2) ) - sqrt( V*delta/(gamma^3) + ( V/(2*gamma^2) )^2 )  
  U <- runif(x,min=0,max=1)
  idx <- which( U < (delta/(delta+gamma*z1)) )
  z2 <- (delta/gamma)^2 /z1[-idx]
  ret <- numeric(x)
  ret[idx] <- z1[idx]
  ret[-idx] <- z2
  
  return(ret)
}


dIG <- function(x,delta,gamma){
	if( delta <= 0 )
    stop("delta must be positive.")
	if( gamma <= 0 )
    stop("gamma must be positive.")

	dens <- delta*exp(delta*gamma)*(x^(-3/2))*exp(-((delta^2)/x+x*gamma^2)/2)/sqrt(2*pi)
	dens
}


## (Multivariate) Normal inverse Gaussian

rNIG <- function(x,alpha,beta,delta,mu,Lambda){
  ## Error check
  #print(match.call())
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
   Lambda <- NA

  if(is.na(Lambda)){
    ## univariate case
    gamma <- sqrt(alpha^2 - beta^2)
    if(gamma <0){
      stop("alpha^2-beta^2 must be positive.")
    }
    
    if (gamma == 0) {
      V = rnorm(x)^2
      Z = delta * delta/V
      X = sqrt(Z) * rnorm(x)
    }else{ 
      Z <- rIG(x,delta,gamma)
      N <- rnorm(x,0,1)
      X <- mu + beta*Z + sqrt(Z)*N
    }
    return(X)
    
  }else{  ## multivariate case
  	
	if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }

    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
    if(tmp <0){
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative.")
    }
    gamma <- sqrt(tmp)
  }
  
  tau <- rIG(x,delta,gamma)
  eta <- rnorm(x*length(beta))
  sqrt.L <- svd(Lambda)
  sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)

  z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
  X <- z
  #  print(str(X))
  return(X)
}


dNIG <- function(x,alpha,beta,delta,mu,Lambda){
  ## Error check
  #print(match.call())
  if(length(alpha)>1||length(delta)>1)
    stop("alpha and delta must be positive reals.")
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if( alpha < 0 )
    stop("alpha must be nonnegative.")
  if( delta <= 0 )
    stop("delta must be positive.")
  if(missing(Lambda))
    Lambda <- NA
  if(is.na(Lambda)){
    #univraiate case
    if(length(beta)>1||length(mu)>1)
      stop("beta and mu must be numeric")
    if( alpha^2 - beta^2 < 0 )
      stop("alpha^2 - beta^2 must be nonnegative.")
    dens <- alpha*delta*exp(delta*sqrt(alpha^{2}-beta^{2})+beta*(x-mu))*besselK(alpha*sqrt(delta^{2}+(x-mu)^{2}),1)/(pi*sqrt(delta^{2}+(x-mu)^{2}))
    dens
  }else{
    #multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    if(nrow(Lambda)!=length(mu)){
      stop("Dimension of Lambda and mu must be equal.")
    }

	if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }
    
    tmp<-alpha-t(beta)%*%Lambda%*%beta
    if(tmp <0){
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative.")
    }
    Lambdainv<-solve(Lambda)
    dens<- alpha^{(nrow(Lambda)+1)/2}*delta*exp(delta*sqrt(tmp)+t(beta)%*%(x-mu))*besselK(alpha*sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)),nrow(Lambda))/(2^{(nrow(Lambda)-1)/2}*pi^{(nrow(Lambda)+1)/2}*(sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)))^{(nrow(Lambda)+1)/2})
    dens
  }	
}

## ## One-dim. NIG
## rNIG <- function(x,alpha=1,beta=0,delta=1,mu=0){
##   gamma <- sqrt(alpha^2 - beta^2)
##   if (gamma == 0) {
##     V = rnorm(x)^2
##     Z = delta * delta/V
##     X = sqrt(Z) * rnorm(x)
##   }else{ 
##     Z <- rIG(x,delta,gamma)
##     N <- rnorm(x,0,1)
##     X <- mu + beta*Z + sqrt(Z)*N
##   }
##   return(X)
## }


## Univariate non-Gaussian stable: corrected Weron's algorithm incorporated
## (20160114, HM) "dstable" still unavailable in YUIMA... incorporate "stabledist" package?

rstable <- function(x,alpha,beta,sigma,gamma){
	if( alpha <= 0 || alpha > 2)
	stop("The index alpha must take values in (0,2].")
	if( beta < -1 || beta > 1)
	stop("The skeweness beta must take values in [-1,1].")
	if( sigma <= 0 )
	stop("The scale sigma must be positive.")
	
  a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
  b <- atan(beta*tan(alpha*pi/2))/alpha
  
  u <- runif(x,-pi/2,pi/2)
  v <- rexp(x,1)
  
  if(alpha!=1){
    s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
	X <- sigma * s + gamma * rep(1,x)
  }else{
    s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log((pi/2)*v*cos(u)/(beta*u + pi/2)))
	X <- sigma * s + (gamma + (2/pi)*beta*sigma*log(sigma)) * rep(1,x)
  }
  
  return(X)
}

## rstable <- function(x,alpha,beta,sigma, gamma, eps){
##   a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
##   b <- atan(beta*tan(alpha*pi/2))/alpha

##   u <- runif(x,-pi/2,pi/2)
##   v <- rexp(x,1)
  
##   if(alpha!=1){
##     s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
##   }else{
##     s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log(v*cos(u)/(beta*u + pi/2)))
##   }
  
##   X <- (eps^(1/alpha)) * sigma * s + gamma * eps * rep(1,x)
##   return(X)
## }



## Positive exponentially tempered stable (AR method)
## ## This must be re-coded later!! (temporarily, rather inefficient)
rpts <- function(x,al,a,b) {
    if( al <= 0 | al>= 1 )
    stop("al must lie in (0,1).")
    if( a <= 0 )
    stop("a must be positive value.")
    if( b <= 0 )
    stop("b must be positive value.")

	sig <- (-a*gamma(-al)*cos(pi*al/2))^(1/al)
	y <- c()
	i <- 1
	while (i <= x) {
		u <- runif(1)
		v <- rstable(1,al,1,sig,0)
		w <- exp(-b*v)
		if (u < w ){
			y <- append(y,v)
			i <- i+1
		}
	}
	return(y)
}


