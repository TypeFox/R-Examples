cal_u <-
function(clsLevel, var_u, var_e, homo, plmm){
  i<-(plmm$clsLev==clsLevel)
  ni<-plmm$ni[i]
  #sum_y_xb_gamma=sum(plmm$y_Xb_gamma0[i][[1]])
  sum_v<-sum(plmm$v0[i][[1]])
  #u_hat=(var_u/(ni*var_u+var_e))*sum_y_xb_gamma
  u_hat<-(var_u/(ni*var_u+var_e))*sum_v
  names(u_hat)<-NULL
  return(u_hat)
}
