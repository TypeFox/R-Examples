FC_xixi <-
function(clsLev, p, plmm){
  i<-(plmm$clsLev==clsLev)   
  ni<-plmm$ni[i]       
  X_tilde<-matrix(plmm$X_tilde0[i][[1]], ncol=p) 
### sum of nixx
  if(p==1){
    sum_nixx<-(ni^2)*mean(X_tilde)^2
  }else{
    x_tilde_bar<-colMeans(X_tilde)
    sum_nixx<-(ni^2)*x_tilde_bar%*%t(x_tilde_bar)
  }
  return(sum_nixx=sum_nixx)
}
