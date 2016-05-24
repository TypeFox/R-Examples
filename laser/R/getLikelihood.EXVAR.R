`getLikelihood.EXVAR` <-
function(mu0, z, lam, SpeciationTimes, Tmax)
{

  x <- SpeciationTimes;
  L2 <- 0; L3 <- 0; 
  L1<- sum(log(2:(length(x) -1))) + (length(x)-2)*log(lam); 
  L2 <- (2*log(pTt.EXVAR(mu0, z, lam, x[1], Tmax)) + 2*rhoFxn.EXVAR(mu0, z, lam, x[1], Tmax));
  for (i in 3:length(x)){
  	L3 <- (L3 + 2*log(pTt.EXVAR(mu0, z, lam, x[i], Tmax)) + rhoFxn.EXVAR(mu0, z, lam, x[i], Tmax));
  }
  #cat(L1, L2, L3, sep='\n')
  Likelihood <- L1 + L2 + L3;
  Likelihood; 

}

