# ' Gibbs sampler for INDEPENDENT Covariates
# '  @param X est la ligne concernee (X[i,]) 
# ' @param components est le vecteur des classes pour i (vecteur creux de taille p avec pr elements non nuls)
# ' @param mixmod is mixmod$details
Gibbs_X_ij_IF<-function(Z=Z,X=X,p=p,mui=mui,sigmai=sigmai,Sigma=Sigma,alpha=alpha,mixmod=mixmod,j=j,components=components,i=i){
   Sigma_j_reste=Sigma[j,-j]
   Sigma_reste_reste=Sigma[-j,-j]
   prodmat=Sigma_j_reste%*%solve(Sigma_reste_reste)
   mu=mui[j]+prodmat%*%(X[-j]-mui[-j])
   sigma=sigmai[j]-prodmat%*%Sigma_j_reste;
   sigmai[j]-Sigma_j_reste%*%solve(Sigma_reste_reste)%*%Sigma_j_reste
   sigma=as.numeric(sigma)
#    print(paste("sigma",sigma))
   if(as.numeric(sigma)<=0){
      print(paste("sigmas<0",sigma,i, j,"try to scale the dataset"));
      sigma=-as.numeric(sigma)
#       stop("bullshit")
   }
   res=rnorm(1,mean=as.numeric(mu),sd=as.numeric(sqrt(sigma)))
#    print(paste("res",res))
   if(is.na(res)){print(sigma);print('NA');res=X[j]}
   return(res)
}