"sde.sim.ea" <-
function(X0, t0, Dt, N, d1, d1.x, k1, k2, phi, max.psi, rh, A){
 
   psi <- function(x) 0.5*d1(1,x)^2 + 0.5*d1.x(1,x)

   if(missing(k1)){
    message("k1 missing, trying numerical minimization...")
    k1 <- optimize(psi, c(0, max.psi))$obj
    message(sprintf("(k1=%5.3f)\n",k1))
   }
   if(missing(k2)){
    message("k2 missing, trying numerical maximization...")
    k2 <- optimize(psi, c(0, max.psi),maximum=TRUE)$obj
    message(sprintf("(k2=%5.3f)\n",k2))
  }
   
  if(missing(phi))
    phi <- function(x) 0.5*d1(1,x) + 0.5*d1.x(1,x) - k1
  else
   phi <- function(x) eval(phi)

  M <- k2-k1
  if(M==0)
   stop("`k1' = `k2' probably due to numerical maximization")

  if(Dt>1/M)
    stop(sprintf("discretization step greater than 1/(k2_k1)"))
   
   if(missing(A))
     A <- function(x) integrate(d1, 0, x)

  if(missing(rh)){
   rh <- function(){
    h <- function(x) exp(A(x) - x^2/(2*Dt))
    f <- function(x) h(x)/dnorm(x,sqrt(Dt))
    maxF <- optimize(f,c(-3*Dt, 3*Dt),maximum=TRUE)$obj
    while(1){
     y <- rnorm(1)
     if( runif(1) < f(y)/maxF )
      return(y)
   }
  }
 }

  x0 <- X0
  X <- numeric(N)
  X[1] <- X0
  rej <- 0
  j <- 1
  while(j <= N){
   y <- x0+rh()
   k <- rpois(1,M*Dt)
   if(k>0){
    t <- runif(k)*Dt
    v <- runif(k)*M
    idx <- order(t)
	t <- c(0, t[idx], Dt)
    v <- v[idx]

	DT <- t[2:(k+2)] - t[1:(k+1)]
	W <- c(0,cumsum(sqrt(DT) * rnorm(k+1))) 
    Y <- x0 + W -(W[k+2] -y+x0)*t/Dt
    if( prod(phi(Y[2:(k+1)]) <= v) == 1){
     j <- j+1
     x0 <- Y[k+2]
     X[j] <-  Y[k+2]
    } else {
     rej <- rej +1
   }
  }
 }
 cat(sprintf("rejection rate: %5.3f\n",rej/(N+rej)))
 X
}

