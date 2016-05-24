# ' Calcul bic local
# ' 
calcul_BIC_theta<-function(theta=theta,X=X,Y=Y,k=k){#ici k est la complexit? totale
  return(-2*log_likelihood(theta=theta,Y=Y,X=X)+(k*log(length(Y))))
}