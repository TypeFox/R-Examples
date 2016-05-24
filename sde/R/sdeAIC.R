sdeAIC <- function(X, theta, b, s, b.x, s.x, s.xx, B, B.x, H, S, guess, ...){
	n <- length(X)
	DELTA <- deltat(X)
	 
	if(missing(theta) || is.null(theta)){
	 if(missing(guess))
	  stop("cannot estimate the model. Specify initial guess values for theta")
	 g1 <- function(theta,X,drift,sigma){
      sum(log(sigma(X[-n],theta) + (diff(X)-DELTA*drift(X[-n],theta))^2/(2*DELTA*sigma(X[-n],theta)^2)))
     }
	 cat("estimating the model...\n")
	 est <- optim(guess,g1,drift=b,sigma=s,X=X,...)
	theta <- est$par
	print(theta)
	}
	
	if(missing(B))
		B <- function(x,theta) b(x,theta)/s(x,theta) - 0.5*s.x(x,theta)
	 if(missing(B.x) && !( missing(b.x) || missing(s.x) || missing(s.xx))){
	  B.x <- function(x,theta) (b.x(x,theta)/s(x,theta) 
	     - b(x,theta)*s.x(x,theta)/(s(x,theta)^2)
	      -0.5*s.xx(x,theta))
	} else {
          stop("error")
	} 
	  
	C1 <- function(x) (B(x,theta)^2)/3 + 0.5*B.x(x,theta)*s(x,theta)
	g <- function(x,y) -0.5*( C1(x) + C1(y) + B(x,theta)*B(y,theta)/3)
	h <- function(x) B(x,theta)/s(x,theta)
	if(missing(H))
	 H <- function(x,y)
		sapply(1:length(x), function(i) integrate(h, x[i], y[i])$value) 
	s1 <- function(x) 1/s(x,theta)
	if(missing(S))
	 S <- function(x,y) 
	  sapply(1:length(x), function(i) integrate(s1, x[i], y[i])$value)
	 
    	u <- function(x,y,theta) (-0.5*log(2*pi*DELTA) - log(s(y,theta)) 
	    - S(x,y)^2/(2*DELTA)  + H(x,y) + DELTA*g(x,y))
	cat("AIC value:\n")
    	-2*sum(u(X[1:(n-1)],X[2:n],theta)) + 2*length(theta)
}
