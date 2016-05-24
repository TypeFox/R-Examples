lk_obs_score <- function(par_comp,lde,lpar,lga,S,R,yv,k,rm,l,J,fv,link,disc,indga,glob,refitem,miss,ltype,
                         XXdis,Xlabel,ZZ0,fort){

# preliminaries
	ncov = dim(XXdis)[2]
	ns = length(Xlabel)
	cov = TRUE
	if(glob) logit_cov = "g" else logit_cov = "m"
# separate parameters
	de = par_comp[1:lde]
	par = par_comp[lde+(1:lpar)]	
	if(disc==1) ga = par_comp[lde+lpar+(1:lga)]			
# Compute log-likelihood
	Piv = prob_multi_glob(XXdis,logit_cov,de,Xlabel)$P
	if(disc==0) ZZ = ZZ0
	if(disc==1){
		if(rm<J){
			gac = rep(1,J); gac[indga] = ga
			ZZ = ZZ0
   			for(j in 1:J){
    			ind = (refitem==j)
	    		ZZ[,,ind] = ZZ[,,ind]*gac[j]
			}
		}
   	}
	P = prob_multi_glob(ZZ,ltype,par)$P
	Phi = array(t(P),c(l,J,k))
	Psi = matrix(1,ns,k)
	if(miss){
		for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
	}else{
		if(fort){
			o = .Fortran("lk_obs",J,as.integer(k),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
           Psi = o$Psi
		}else{
   		    for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
   		}            	
	}
	if(k==1) Pj=Psi else Pj = Psi*Piv
	pm = rowSums(Pj)
	lk = sum(yv*log(pm))
#	print(lk)
# ---- E-step ----
	V = ((yv/pm)%o%rep(1,k))*Piv*Psi; sV = colSums(V)
# ---- M-step ----
	YY = matrix(0,J*k,l)
	count = 0
	for(c in 1:k) for(j in 1:J){
		count = count+1
		for(y in 1:l){
			ind = (S[,j]==(y-1))
			if(miss) YY[count,y] = sum(V[ind,c]*R[ind,j]) else YY[count,y] = sum(V[ind,c])			
		}
	}
	if(disc==0){
		sc_ga = NULL
	}else{
		if(rm<J){
			ZZ1 = array(0,c(l-1,J,J*k))
			count = 0
			for(c in 1:k) for(j in 1:J){
				count = count+1
				ZZ1[,j,count] = ZZ0[,,count]%*%par
			}
			dimz = dim(ZZ1)
			dimz[2] = dimz[2]-length(fv)
			if(rm==1){
				ZZ1int = ZZ1[,fv,]	
			}else{
				ZZ1int = apply(ZZ1[,fv,],c(1,3),sum)							
			}
			ZZ1 = array(ZZ1[,-fv,],dimz)
			if(l==2) ZZ1int = matrix(ZZ1int,1,length(ZZ1int))				
			sc_ga = est_multi_glob(YY,ZZ1,ltype,be=ga,Int=ZZ1int,only_sc=TRUE)$sc
		}
		ZZ = ZZ0
    	for(j in 1:J){
			ind = (refitem==j)
	    	ZZ[,,ind] = ZZ[,,ind]*gac[j]
		}
	}
	sc_par = est_multi_glob(YY,ZZ,ltype,be=par,only_sc=TRUE)$sc
# Update piv
	sc_de = est_multi_glob(V,XXdis,logit_cov,Xlabel,de,only_sc=TRUE)$sc
# output
	sc = c(sc_de,sc_par,sc_ga)
	out = list(lk=lk,sc=sc)
	
}