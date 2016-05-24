GCV <-
function(h, yXb, d, N, poly.index, plmm, sd_T){
  h<-exp(h)
  if(d==2){h<-h*sd_T}  
  gamma_S<-sapply(as.list(1:N), .gamma_S, yXb=yXb, d=d, h=h, N=N, poly.index=poly.index, plmm=plmm)

  gamma<-unlist(gamma_S[rownames(gamma_S)=="gamma.k"])
#gamma_sm=sm.regression(y=yXb, x=plmm$T_mat, h=h, eval.grid=F, eval.points=plmm$T_mat, poly.index=poly.index, display="none")$estimate

  trS<-sum(unlist(gamma_S[rownames(gamma_S)=="S.kk"]))
  GCV<-sum((yXb-gamma)*(yXb-gamma))/(N*(1-trS/N)^2)  
  return(c(GCV))                     
}
