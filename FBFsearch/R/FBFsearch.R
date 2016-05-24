FBF_LS=function(Corr,nobs,G_base,h,C,n_tot_mod) {
.Call("FBF_LS", as.matrix(Corr),as.vector(nobs),as.matrix(G_base),as.vector(h),as.vector(C),as.vector(n_tot_mod))
}

FBF_RS=function(Corr,nobs,G_base,h,C,n_tot_mod,n_hpp) {
.Call("FBF_RS", as.matrix(Corr),as.vector(nobs),as.matrix(G_base),as.vector(h),as.vector(C),as.vector(n_tot_mod),as.vector(n_hpp))
}

FBF_GS=function(Corr,nobs,G_base,h,C,n_tot_mod,n_hpp) {
.Call("FBF_GS", as.matrix(Corr),as.vector(nobs),as.matrix(G_base),as.vector(h),as.vector(C),as.vector(n_tot_mod),as.vector(n_hpp))
}

