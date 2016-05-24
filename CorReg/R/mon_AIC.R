mon_AIC<-function(theta=theta,Y=Y,X=X,intercept=TRUE){
   #theta est le vecteur A sans l'ecart-type qu'on calcule donc ici
   if(intercept){X=cbind(1,X)}
   sigma=sd(Y-X%*%theta)
   theta=c(sigma,theta)
   AIC=-2*log_likelihood(theta=theta,Y=Y,X=X)+2*length(theta[theta!=0])
   return(AIC)
}