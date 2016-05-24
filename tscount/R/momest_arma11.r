momest_arma11 <- function(ts){
  #Moment estimators for ARMA(1,1) process
  #ts: Numeric vector. Time series assumed to be a relaisation of an ARMA(1,1) process.
  m <- mean(ts)
  a <- acf(ts, lag.max=2, plot=FALSE)$acf[-1,,1]
  psi1 <- a[2]/a[1]
  #The MA parameter is the solution of a quadratic equation in monic form (i.e. the first coefficient is 1) which is solved using the pq formula:
  Qe <- (psi1^2-2*a[1]*psi1+1)/(a[1]-psi1)
  Pu <- 1
  discriminant <- Qe^2-4*Pu
  if(discriminant>=0){ #the solutions are real
    solutions <- (-Qe + c(+1,-1)*sqrt(discriminant))/2
  }else{ #the soulutions are complex (namely -Qe/2 + c(+1,-1)*i*sqrt(-discriminant)), use their orthogonal projection on the real axis instead
    solutions <- -Qe/2
  }
  theta1 <- solutions[abs(solutions)<=1][1] #choose a solution for which the resulting ARMA(1,1) process can be stationary (although this is not yet guaranteed only by fulfilling abs(.)<=1)
  sigmasq <- (1-psi1^2)/(1-theta1^2)*var(ts)
  result <- c(ar1=psi1, ma1=-theta1, intercept=m) #different parametrisation
  return(result)
}
