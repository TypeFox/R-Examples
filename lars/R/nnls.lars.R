nnls.lars <-
function(active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE, use.Gram = TRUE)
{
### Modified 05/15/03 to allow for more than one addition to the set
### Lawson and Hanson page 161
### Go back to the first positive coefficent vector; can assume its in order
### Note that X'y is constant for all these guys;
### we assume WOLOG this constant is 1
### We also assume we have come into this because we have a negative coeff  
### If use.Gram is FALSE, then Gram comes in as x
  if(!use.Gram) x <- Gram	# to avoid confusion
  M<-m <- length(active)
  im <- seq(m)
  positive <- im
  zero <- NULL
  ### Get to the stage where beta.old is all positive
  while(m>1) {
    zero.old<-c(m,zero)
    R.old <- downdateR(R, m)
    beta0 <- backsolve(R.old, backsolvet(R.old, Sign[ - zero.old]))*Sign[-zero.old]
    beta.old <- c(beta0,rep(0,length(zero.old)))
    if(all(beta0 >0))break
    m <-m-1
    zero<-zero.old
    positive<-im[-zero]
    R<-R.old
    beta<-beta.old
  }
### Now we do the NNLS backtrack dance
  while(TRUE) {
    while(!all(beta[positive] > 0)) {
      alpha0 <- beta.old/(beta.old - beta)
      alpha <- min(alpha0[positive][(beta <= 0)[positive]])
      beta.old <- beta.old + alpha * (beta - beta.old)
      dropouts<-match(alpha,alpha0[positive],0)
### Used to have the following line, but this failed occasionally
###   dropouts <- seq(positive)[abs(beta.old[positive]) < eps]
      for(i in rev(dropouts)) R <- downdateR(R, i)
      positive <- positive[ - dropouts]	
                                        # there is an order in R
      zero <- im[ - positive]
      beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * 
        Sign[positive]
      beta <- beta.old * 0
      beta[positive] <- beta0
    }
### Now all those in have a positive coefficient
    if(use.Gram) {
      w <- 1 - Sign * drop(Gram %*% (Sign * beta))	
                                        #should be zero for some
    }
    else {
      jw <- x %*% (Sign * beta)
      w <- 1 - Sign * drop(t(jw) %*% x)
    }
    if((length(zero) == 0) || all(w[zero] <= 0))
      break
    add <- order(w)[M]
    if(use.Gram) {
      R <- updateR(Gram[add, add], R, drop(Gram[add, 
                                                positive]), Gram = TRUE,eps=eps)
    }
    else {
      R <- updateR(x[, add], R, x[, positive], Gram = FALSE,eps=eps)
    }
    positive <- c(positive, add)
    zero <- setdiff(zero, add)
    beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * Sign[
                                                        positive]
    beta[positive] <- beta0
  }
  if(trace)
    {
      dropouts<-active[-positive]
      for(i in dropouts){
          cat("NNLS Step:\t Variable", i, "\tdropped\n")
        }
    }
  list(active = active[positive], R = R, beta = beta, positive = positive
       )
}

