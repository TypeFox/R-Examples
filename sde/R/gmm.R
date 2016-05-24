gmm <- function(X, u, dim, guess, lower, upper,maxiter=30,tol1=1e-3,tol2=1e-3){
  if(!is.ts(X)) 
	stop("Please provide a `ts' object")
  DELTA = deltat(X)
  n <- length(X) 

  if(!missing(guess)){
   if(length(guess)>0){
    dim <- length(guess)
    cat(sprintf("\nDimension of parameter space set to %d", dim))
   }
  }
   
  if(missing(dim))
   stop("Please specify dimension of parameter space")  

  H <-function(theta)
   apply(u(X[2:n], X[1:(n-1)], theta, DELTA), 2, mean)

  Q <-function(theta) sum(H(theta)^2)

  S <- function(j, theta)
   ( (t(u(X[(j+2):n],X[(j+1):(n-1)], theta, DELTA)) %*%  
     u(X[2:(n-j)],X[1:(n-j-1)], theta, DELTA))/n )
 
  ell <- n-2
  w <- 1-(1:ell)/(ell+1) # Bartlet weights

  cat("\nInitial values for the optimization algorithm ")
  if(missing(guess)){ 
    guess <- runif(dim)
	cat("(random)\n")
  } else {
   cat("\n")
  }
  print(guess)
  
  if(missing(lower)) 
	lower <- rep(-Inf, dim)
  if(missing(upper)) 
	upper <- rep(Inf, dim)

  cat("\nOptimization contraints\n")
  ct <- as.matrix(cbind(lower, upper))
    rownames(ct) <- paste("theta",1:dim,sep="")
    colnames(ct) <- c("lower", "upper")
    print(ct)
    cat("\nRunning optimizer...\n")


  theta1 <- optim(guess,Q,method="L-BFGS-B", upper=upper, lower=lower)$par
  for(i in 1:dim){
   names(theta1)[i] <- sprintf("theta%d",i)
  } 
  cat("\nFirst stage estimates:\n")
  print(theta1)
  
  cat("\nStarting second stage...\n") 
  goOn <- TRUE
  iter <- 0
  while(goOn){
   iter <- iter + 1
   S.hat <- S(0, theta1)
   for(i in 1:ell)
    S.hat = S.hat + w[i]*(S(i,theta1)+t(S(i,theta1)))
   W <- solve(S.hat)
   Q1 <-function(theta) H(theta) %*% W %*% H(theta)
   fit <- optim(theta1,Q1,method="L-BFGS-B", upper=upper, lower=lower, hessian=TRUE)
   theta2 <- fit$par
   val <- fit$value
   hes <- fit$hessian
   out <- c(theta2, val, sum(abs(theta1-theta2)))
   names(out) <- c(names(theta1),"Q1","|theta1-theta2|")
   print(out)
   names(theta2) <- names(theta1)  
   if(sum(abs(theta1-theta2))<tol1 || val<tol2 || iter>maxiter)
    goOn <- FALSE
   theta1 <- theta2
 } 
 
  list(par=theta1, val=val, hessian=hes)
}
