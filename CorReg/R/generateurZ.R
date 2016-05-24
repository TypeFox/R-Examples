# ' generateur de Z aleatoire
# MH seulement mais all?g? pour ?tre plus rapide
generateurZ<-function(p=p,maxit=NULL){
  if(is.null(maxit)){
    #maxit=round(p^2/4) 
    maxit=floor(runif(1,min=0,max=p^2))
  }
  Z=Matrix(0,ncol=p,nrow=p)
  candidats=(1:p^2)[-seq(1,p^2,p+1)] #tout sauf la diagonale
  #on initialise
  for(i in 1:maxit){
    qui=arrayInd(sample(candidats,size=1),c(p,p))#ordre al?atoire des candidats
    i=qui[1]
    j=qui[2]
    Z[,i]=Z[i,j]*Z[,i]#annule si ? 0 (donc si on le passe ? 1) neutre sinon
    Z[j,]=Z[i,j]*Z[j,]#annule si ? 0 (donc si on le passe ? 1) neutre sinon
    Z[i,j]=1-Z[i,j]#on modifie le point candidat
  }
  return(Z)
}