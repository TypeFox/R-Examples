polarisation.ER <-
function(alpha,rho,comp=FALSE){
  means<-rho[,1]
  if(names(rho)[1]!="means") print("Warning: First column of rho is not called means. Ensure that it entails are the groups' mean incomes!")
  if(names(rho)[2]!="shares") print("Warning: First column of rho is not called shares. Ensure that it entails are the groups' population shares!")
  shares    <- rho[,2]
  sharesi   <- shares^(1+alpha)
  sharesj   <- shares
  meansi    <- means%*%t(rep(1,length(means)))
  meansj    <- t(means%*%t(rep(1,length(means))))
  meansdiff <- abs(meansi-meansj)
  ER    <- sharesi%*%meansdiff%*%sharesj
  if(comp){
    ERcomp<-matrix(NA,length(means),length(means))
    for(i in 1:length(means)){
      for(j in 1:length(means)){
        ERcomp[i,j] <- rho[i,2]^(1+alpha)*rho[j,2]*abs(rho[i,1]-rho[j,1])
      }
    }
    list(P=ER,ERcomp=ERcomp,means=means,shares=shares) 
  }else list(P=ER,means=means,shares=shares)   
}
