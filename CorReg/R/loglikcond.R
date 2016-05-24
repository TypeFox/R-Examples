
# ' X is the entire dataset (const &)

loglikcond<-function(X=X,mui=mui,Sigma=Sigma,M=M,i=i,Zc=Zc){ 
   res=dmvnorm(x=X[i,],mean=mui,sigma=Sigma)#vraisemblance globale de l'individu
   quimank_i=which(M[i,]==1)
   if(length(quimank_i)>0){
      Sigma_M_obs=matrix(Sigma[quimank_i,-quimank_i],nrow=length(quimank_i))
      Sigma_obs_obs=Sigma[-quimank_i,-quimank_i]
      prodmat=Sigma_M_obs%*%solve(Sigma_obs_obs)
      mu=mui[quimank_i]+prodmat%*%(X[i,-quimank_i]-mui[-quimank_i])
      sigma=Sigma[quimank_i,quimank_i]-prodmat%*%t(Sigma_M_obs)#matrice variance-covariance
      res=res/dmvnorm(x=X[i,quimank_i],mean = mu,sigma=sigma)#on divise par la vraisemblance conditionnelle des manquants pour obtenir la marginale
      #attention sigma minuscule dans les deux lignes precedentes  
   }
   return(log(res))
}

