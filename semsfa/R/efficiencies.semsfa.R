efficiencies.semsfa<-function(semobj,log.output=TRUE,...){
    z_sem=(semobj$y-semobj$fitted)*semobj$lambda/semobj$sigma
    if(semobj$ineffDecrease)
      {
      u_sem=round((semobj$sigma*semobj$lambda)/(1+semobj$lambda^2) *(dnorm(z_sem)/(1-pnorm(z_sem)) - z_sem),4)
      if(log.output){efficiencies=exp(-u_sem)} else {efficiencies=(semobj$fitted-u_sem)/semobj$fitted}
      }
    else
      {
      u_sem=round((semobj$sigma*semobj$lambda)/(1+semobj$lambda^2) *(dnorm(z_sem)/(pnorm(z_sem)) + z_sem),4)
      if(log.output){efficiencies=exp(-u_sem)} else {efficiencies=(u_sem-semobj$fitted)/semobj$fitted}
      }
    semobj$u<- u_sem
    semobj$efficiencies<-efficiencies
    return(semobj)
}