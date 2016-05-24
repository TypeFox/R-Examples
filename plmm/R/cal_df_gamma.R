cal_df_gamma <-
function(N, d, h, plmm, poly.index){ 
  S_SS<-sapply(as.list(1:N), cal_S_SS, d=d, h=h, N=N, poly.index=poly.index, plmm=plmm)
  
  trS<-sum(unlist(S_SS[rownames(S_SS)=="S.kk"]))
  trSS<-sum(unlist(S_SS[rownames(S_SS)=="S.kj2"]))
  
  return(c(2*trS-trSS))
}
