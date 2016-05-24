cal_S_SS <-
function(k, d, h, N, plmm, poly.index){
  T_star<-plmm$T_mat-matrix(1, nrow=N, ncol=d)%*%diag(plmm$T_mat[k,], ncol=d)
  K_multi<-exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2)) 
  K<-apply(K_multi, 1, prod)
  
  if(poly.index==0){ # NW
    T_star<-matrix(1, nrow=N, ncol=1)
  }else{T_star<-cbind(1, T_star)} #Local Linear
  
  KT<-T_star*K
  TKT<-t(T_star)%*%KT
  TKTinv<-solve(TKT)
  s.k<-(TKTinv%*%t(KT))[1,] # k-th row of S

  S.kk<-s.k[k]# [S]_kk
  
### tr(SS)=sum{([S]_kj)^2}
  S.kj2<-sum(s.k^2) # sum of squares of the k-th row of S
  
  return(list(S.kk=S.kk, S.kj2=S.kj2))
}
