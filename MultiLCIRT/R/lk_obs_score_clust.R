lk_obs_score_clust <- function(par_comp,lde1,lde2,lpar,lga,S,R,kU,kV,rm,l,J,fv,link,disc,indga,refitem,miss,ltype,
                               WWdis,Wlabel,XXdis,Xlabel,ZZ0,clust,fort){

# preliminaries
	ncov = dim(XXdis)[2]
	ns = nrow(S)
	cov = TRUE
	nclust = max(clust)	
# separate parameters
	de1 = par_comp[1:lde1]
	de2 = par_comp[lde1+(1:lde2)]
	par = par_comp[lde1+lde2+(1:lpar)]	
	if(disc==1) ga = par_comp[lde1+lde2+lpar+(1:lga)]			
# Compute log-likelihood
	if(kU==1){
		La = matrix(1,nclust,1)
	}else{
	    out = prob_multi_glob(WWdis,"m",de1,Wlabel)
	    Ladis = out$Pdis; La = out$P
	}
	if(kV==1){
		Piv = rep(1,ns)
	}else{
	    out = prob_multi_glob(XXdis,"m",de2,Xlabel)
    	Pdis = out$Pdis; Piv = out$P
    }
 	Piv = array(t(Piv),c(kV,ns,kU))
 	Piv = aperm(Piv,c(2,1,3))
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
	Phi = array(t(P),c(l,J,kV))
	Psi = matrix(1,ns,kV)
	if(miss){
		for(j in 1:J) for(c in 1:kV) Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
	}else{
		if(fort){
			o = .Fortran("lk_obs",J,as.integer(kV),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
           Psi = o$Psi
		}else{
   		    for(j in 1:J) for(c in 1:kV) Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
   		}            	
	}
	if(kV==1){
		Pj = pm = Psi; lk = sum(log(pm))
	}else{
		lPj = matrix(0,nclust,kU)
		for(u in 1:kU) for(cl in 1:nclust){
			ind = which(clust==cl)
			if(length(ind)==1){
				lPj[cl,u] = log(Psi[ind,]%*%Piv[ind,,u])+log(La[cl,u])
			}else{
				lPj[cl,u] = sum(log(rowSums(Psi[ind,]*Piv[ind,,u])))+log(La[cl,u])
			}
		}
		mlPj = apply(lPj,1,max)
		pm = rowSums(exp(lPj-mlPj))
		lk = sum(log(pm)+mlPj)
	}
#	print(lk)
# ---- E-step ----
	if(kV==1){
		Vclust = matrix(1,nclust,1)
		V = matrix(1,ns,1)
	}else{
		Tmp = exp(lPj-mlPj)
		Vclust = Tmp*(1/rowSums(Tmp))
		Vcomp = NULL; V = 0
		Tmp1 = NULL
		for(u in 1:kU){
			Tmp = Psi*Piv[,,u]
			Tmp = Tmp*(Vclust[clust,u]/rowSums(Tmp))
			V = V+Tmp
			Vcomp = rbind(Vcomp,Tmp)
			Tmp1 = rbind(Tmp1,colSums(Tmp)/sum(Tmp))
		}
	}
	sV = colSums(V)
# ---- M-step ----
	YY = matrix(0,J*kV,l)
	count = 0
	for(c in 1:kV) for(j in 1:J){
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
			ZZ1 = array(0,c(l-1,J,J*kV))
			count = 0
			for(c in 1:kV) for(j in 1:J){
				count = count+1
				ZZ1[,j,count] = ZZ0[,,count]%*%par
			}
			dimz = dim(ZZ1)
			dimz[2] = dimz[2]-length(fv)
			if(rm==1){
				ZZ1int = ZZ1[,fv,]	
			}else{
				ZZ1int = apply(array(ZZ1[,fv,],c(l-1,rm,J*kV)),c(1,3),sum)
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
	if(kU==1) sc_de1 = NULL
	else sc_de1 = est_multi_glob(Vclust,WWdis,"m",Wlabel,de1,only_sc=TRUE)$sc
	if(kV==1) sc_de2 = NULL
	else sc_de2 = est_multi_glob(Vcomp,XXdis,"m",Xlabel,de2,only_sc=TRUE)$sc	
# output
	sc = c(sc_de1,sc_de2,sc_par,sc_ga)
	out = list(lk=lk,sc=sc)
	
}