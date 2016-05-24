FC_hetero <-
function(N, p, m, plmm, invarIND){
  XQXXQy_axx<-sapply(as.list(plmm$clsLev), cal_XQXXQy_axx, p=p, plmm=plmm, invarIND=invarIND)

  XQy0<-XQXXQy_axx[rownames(XQXXQy_axx)=="XQy"]# List of m (p by 1) vectors
  XQX0<-XQXXQy_axx[rownames(XQXXQy_axx)=="XQX"]# List of m (p by p) matrices
  axx0<-XQXXQy_axx[rownames(XQXXQy_axx)=="axx"]
#browser()

  XQy1<-matrix(unlist(XQy0), nrow=p-sum(invarIND))
  XQy<-rowSums(XQy1)
  XQX1<-matrix(unlist(XQX0), nrow=(p-sum(invarIND))*(p-sum(invarIND)))
  XQX<-matrix(rowSums(XQX1), ncol=p-sum(invarIND))

#browser()  
  axx1<-matrix(unlist(axx0), nrow=p*p)
  axxInv<-solve( matrix(rowSums(axx1), ncol=p) )

  beta_w<-solve(XQX)%*%XQy

### eQe and nixx ###
  eQe_eta_Cfc<-sapply(as.list(plmm$clsLev), cal_eQe_eta_Cfc, beta_w=beta_w, axxInv=axxInv, p=p, plmm=plmm, invarIND=invarIND)
  
## for var_e
  eQe<-eQe_eta_Cfc[rownames(eQe_eta_Cfc)=="eQe"]# List of m scalars
  eQe<-sum(unlist(eQe))

## for var_u
  eta_Cfc<-eQe_eta_Cfc[rownames(eQe_eta_Cfc)=="eta_Cfc"]
  eta_Cfc<-sum(unlist(eta_Cfc))
  
  return(list(eQe=eQe, eta_Cfc=eta_Cfc))
}
