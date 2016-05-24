bootstrap_lm_basic <-
function(piv,Pi,Psi,n,B=100,start=0,mod=0,tol=10^-6){

#        [lk,piv,Pi,Psi,np,aic,bic,lkv] = draw_lm_basic(piv,Pi,Psi,n)
#
# Bootstrap EM estimates piv, Pi and Psi

# Preliminaries
	k = length(piv)
	c = dim(Psi)[1]
	T = dim(Pi)[3]
# Reparametrize
    mPsi = 0; mpiv = 0; mPi = 0
    m2Psi = 0; m2piv = 0; m2Pi = 0
#    mth = 0; m2th = 0;
    for(b in 1:B){
	    	out = draw_lm_basic(piv,Pi,Psi,n)
	    	Sb = out$S; yvb = out$yv
		ns = dim(Sb)[1]
	    	out = est_lm_basic(Sb,yvb,k,start,mod,tol)
	    	mPsi = mPsi+out$Psi/B; mpiv = mpiv+out$piv/B; mPi = mPi+out$Pi/B
	    	m2Psi = m2Psi+out$Psi^2/B; m2piv = m2piv+out$piv^2/B; m2Pi = m2Pi+out$Pi^2/B
	    	mth = mth+out$th/B; m2th = m2th+out$th/B
    }
    sePsi = sqrt(m2Psi-mPsi^2); sepiv = sqrt(m2piv-mpiv^2); sePi = sqrt(m2Pi-mPi^2)
#    seth = sqrt(m2th-mth^2)
    out = list(mPsi=mPsi,mpiv=mpiv,mPi=mPi,sePsi=sePsi,sepiv=sepiv,sePi=sePi)

}
