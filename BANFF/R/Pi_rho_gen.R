#####Generating selecitons of possible pairs of pi and rho
Pi_rho_gen=function(piall,rhoall)
{ 
  pirhopair=list()
  for (j in 1:length(piall))
  {for (k in (j+1):length(rhoall))
  {
    pirhopair$rho0=c(pirhopair$rho0,rep(rhoall[j],length(piall)))
    pirhopair$rho1=c(pirhopair$rho1,rep(rhoall[k],length(piall)))
  }
  }
  for (i in 1:(length(rhoall)*(length(rhoall)-1)/2))
  {pirhopair$pi0=c(pirhopair$pi0,piall[1:length(piall)])
  }
  return(pirhopair)
}
