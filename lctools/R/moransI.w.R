moransI.w <- function(x, w){

    Obs <- length(x)
    n <- Obs
    mean.x <- mean(x)
  
    Xi <- x-mean.x
    Xi.sqr <- Xi*Xi
    Z <- matrix(rep((Xi),Obs),Obs,Obs)
    Z.T <- t(Z)
  
    diag(w)<-0
  
    moran.nom.x <- sum(w*Z*Z.T)
    moran.denom.x <- sum(Xi.sqr)
  
    moran <- (Obs/sum(w))*(moran.nom.x/moran.denom.x)
  
  
    #Expected I E(I)
    E.I <- (-1)/(Obs-1)
  
    S0 <- sum(w)
    
    S1 <- (1/2.0) * sum((w + t(w))^2)
    
    S2 <- sum((apply(w, 1, sum) + apply(w, 2, sum))^2)
    
    b2 <- (sum((x - mean.x)^4)/n)/((sum((x - mean.x)^2)/n)^2)
    
    #Var(I)
    Var.I.resampling <- (((n^2) * S1 - n*S2 + 3 * (S0^2))/(((n^2)-1)*(S0^2)))-(E.I^2)
    Var.I.randomization <- (n*((n^2-3*n+3)*S1-n*S2+3*S0^2))/((n-1)*(n-2)*(n-3)*S0^2)-(b2*((n^2-n)*S1-2*n*S2+6*S0^2))/((n-1)*(n-2)*(n-3)*S0^2)-(E.I^2)
    
    Z.I.resampling <- (moran-E.I)/sqrt(Var.I.resampling)
    Z.I.randomization <- (moran-E.I)/sqrt(Var.I.randomization)
    
    pv.resampling <- 2*pnorm(-abs(Z.I.resampling))
    pv.randomization <- 2*pnorm(-abs(Z.I.randomization))
    
    Results <- list(Morans.I=moran, Expected.I=E.I, z.resampling=Z.I.resampling,z.randomization=Z.I.randomization, 
                    p.value.resampling=pv.resampling,p.value.randomization=pv.randomization)
    
}