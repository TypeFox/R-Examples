SA_hetero <-
function(N, p, m, plmm, invarIND){  
  XQXXQy_XPXXPy<-sapply(as.list(plmm$clsLev), cal_XQXXQy_XPXXPy, p=p, plmm=plmm, invarIND)

### Within Beta
  XQy0<-XQXXQy_XPXXPy[rownames(XQXXQy_XPXXPy)=="XQy"]# List of m (p by 1) vectors
  XQy1<-matrix(unlist(XQy0), nrow=p-sum(invarIND))
  XQy<-rowSums(XQy1)
  XQX0<-XQXXQy_XPXXPy[rownames(XQXXQy_XPXXPy)=="XQX"]# List of m (p by p) matrices 
  XQX1<-matrix(unlist(XQX0), nrow=(p-sum(invarIND))*(p-sum(invarIND)))
  XQX<-matrix(rowSums(XQX1), ncol=(p-sum(invarIND)))
  beta_w<-solve(XQX)%*%XQy
  
### Between Beta
  XPy0<-XQXXQy_XPXXPy[rownames(XQXXQy_XPXXPy)=="XPy"]# List of m (p by 1) vectors
  XPy1<-matrix(unlist(XPy0), nrow=p)
  XPy<-rowSums(XPy1)
  XPX0<-XQXXQy_XPXXPy[rownames(XQXXQy_XPXXPy)=="XPX"]# List of m (p by p) matrices 
  XPX1<-matrix(unlist(XPX0), nrow=p*p)
  invXPX<-solve(matrix(rowSums(XPX1), ncol=p))
  beta_b<-invXPX%*%XPy

### Quadratic estimation  ###
  eQe_vPv_eta_Csa<-sapply(as.list(plmm$clsLev), cal_eQe_vPv_eta_Csa, beta_w=beta_w, beta_b=beta_b, invXPX, p=p, plmm=plmm, invarIND=invarIND)
  
## for var_e
  eQe<-eQe_vPv_eta_Csa[rownames(eQe_vPv_eta_Csa)=="eQe"]# List of m scalars
  eQe<-sum(unlist(eQe))

### for var_u
  vPv<-sum(unlist(eQe_vPv_eta_Csa[rownames(eQe_vPv_eta_Csa)=="vPv"]))
  eta_Csa<-sum(unlist(eQe_vPv_eta_Csa[rownames(eQe_vPv_eta_Csa)=="eta_Csa"]))

  return(list(eQe=eQe, vPv=vPv, eta_Csa=eta_Csa))
}
