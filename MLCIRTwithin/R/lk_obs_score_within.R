lk_obs_score_within <- function(par_comp,lde1,lde2,lpar,lga1,lga2,S,R,yv,k1,k2,rm1,rm2,lv,J,fv,link,disc,indga1,indga2,
                                glob,refitem,miss,ltype,XX1dis,XX2dis,Xlabel,ZZ0,fort){

# preliminaries
	lm = max(lv)
	ncov1 = dim(XX1dis)[2]
	ncov2 = dim(XX2dis)[2]
	ns = length(Xlabel)
	cov = TRUE
	if(glob) logit_cov = "g" else logit_cov = "m"
	Aggr1 = diag(k1)%x%matrix(1,1,k2)
	Aggr2 = matrix(1,1,k1)%x%diag(k2)	
	rm = rm1*rm2
# separate parameters
	de1 = par_comp[1:lde1]
	de2 = par_comp[lde1+1:lde2]
	par = par_comp[lde1+lde2+(1:lpar)]	
	if(disc==1){
		ga1 = par_comp[lde1+lde2+lpar+1:lga1]			
		ga2 = par_comp[lde1+lde2+lpar+lga1+1:lga2]
	}			
# Compute log-likelihood
	if(k1>1) Piv1 = prob_multi_glob(XX1dis,logit_cov,de1,Xlabel)$P
	if(k2>1) Piv2 = prob_multi_glob(XX2dis,logit_cov,de2,Xlabel)$P
	if(k1==1) Piv = Piv2
	if(k2==1) Piv = Piv1
    if(k1>1 & k2>1){
   		Piv = matrix(0,ns,k1*k2)
		for(i in 1:ns) Piv[i,] = Piv1[i,]%x%Piv2[i,]
	}
	if(disc==0) ZZ = ZZ0
	if(disc==1){
		ga1c = rep(1,J); ga1c[indga1] = ga1
		ga2c = rep(1,J); ga2c[indga2] = ga2
		ZZ = ZZ0
		for(j in 1:J){
			ind = (refitem==j)
   			ZZ[,1:(k1*rm1),ind] = ZZ[,1:(k1*rm1),ind]*ga1c[j]
   			ZZ[,k1*rm1+1:(k2*rm2),ind] = ZZ[,k1*rm1+1:(k2*rm2),ind]*ga2c[j]
		}
	}
	P = prob_multi_glob_gen(ZZ,ltype,par)$P
	Phi = array(t(P),c(lm,J,k1*k2))
	Psi = matrix(1,ns,k1*k2)
	if(miss){
		for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
	}else{
		# if(fort){
			# o = .Fortran("lk_obs",J,as.integer(k1*k2),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
           # Psi = o$Psi
		# }else{
	    for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
   		# }            	
	}
	if(k1==1 & k2==1) Pj=Psi else Pj = Psi*Piv
	pm = rowSums(Pj)
	lk = sum(yv*log(pm))
#	print(lk)
# ---- E-step ----
	V = ((yv/pm)%o%rep(1,k1*k2))*Piv*Psi; sV = colSums(V)
# ---- M-step ----
	YY = matrix(NA,J*k1*k2,lm)
	count = 0
	for(c in 1:(k1*k2)) for(j in 1:J){
		count = count+1
		for(y in 1:lv[j]){
			ind = (S[,j]==(y-1))
			if(miss) YY[count,y] = sum(V[ind,c]*R[ind,j]) else YY[count,y] = sum(V[ind,c])			
		}
	}
	if(disc==0){
		sc_ga1 = NULL; sc_ga2 = NULL
	}else{
		if(rm<J){
			ZZ1 = ZZ2 = array(NA,c(lm-1,J,J*k1*k2))
			count = 0
			for(c1 in 1:k1) for(c2 in 1:k2) for(j in 1:J){
				count = count+1
				ZZ1[1:(lv[j]-1),,count] = 0
				ZZ1[1:(lv[j]-1),j,count] = ZZ0[1:(lv[j]-1),1:(k1*rm1),count]%*%par[1:(k1*rm1)]
				ZZ2[1:(lv[j]-1),,count] = 0
				ZZ2[1:(lv[j]-1),j,count] = ZZ0[1:(lv[j]-1),k1*rm1+1:(k2*rm2),count]%*%par[k1*rm1+1:(k2*rm2)]
			}
			ZZ1 = array(ZZ1[,indga1,],c(lm-1,length(ga1),J*k1*k2))
			ZZ2 = array(ZZ2[,indga2,],c(lm-1,length(ga2),J*k1*k2))
			ind = (k1*rm1+k2*rm2+1):dim(ZZ0)[2]
			ZZ1Int = ZZ2Int = array(NA,c(lm-1,J*k1*k2))
			count = 0
			for(c1 in 1:k1) for(c2 in 1:k2) for(j in 1:J){
				count = count+1
				ZZ1Int[1:(lv[j]-1),count] = ga2c[j]*ZZ0[1:(lv[j]-1),k1*rm1+1:(k2*rm2),count]%*%par[k1*rm1+1:(k2*rm2)]+ZZ0[1:(lv[j]-1),ind,count]%*%par[ind]
				ZZ2Int[1:(lv[j]-1),count] = ga1c[j]*ZZ0[1:(lv[j]-1),1:(k1*rm1),count]%*%par[1:(k1*rm1)]+ZZ0[1:(lv[j]-1),ind,count]%*%par[ind]
			}
			sc_ga1 = est_multi_glob_gen(YY,ZZ1,ltype,be=ga1,Int=ZZ1Int,only_sc=TRUE)$sc
			sc_ga2 = est_multi_glob_gen(YY,ZZ2,ltype,be=ga2,Int=ZZ2Int,only_sc=TRUE)$sc
		}
		ga1c = rep(1,J); ga1c[indga1] = ga1
		ga2c = rep(1,J); ga2c[indga2] = ga2
		ZZ = ZZ0
		for(j in 1:J){
	    	ind = (refitem==j)
			ZZ[,1:(k1*rm1),ind] = ZZ[,1:(k1*rm1),ind]*ga1c[j]
			ZZ[,k1*rm1+1:(k2*rm2),ind] = ZZ[,k1*rm1+1:(k2*rm2),ind]*ga2c[j]
		}
	}
	# out = est_multi_glob(YY,ZZ,ltype,be=par,Dis=Dis)						
    # par = out$be; P = out$P
	# Phi = array(t(P),c(l,J,k1*k2))
# #######
		# if(rm<J){
			# ZZ1 = array(0,c(l-1,J,J*k))
			# count = 0
			# for(c in 1:k) for(j in 1:J){
				# count = count+1
				# ZZ1[,j,count] = ZZ0[,,count]%*%par
			# }
			# dimz = dim(ZZ1)
			# dimz[2] = dimz[2]-length(fv)
			# if(rm==1){
				# ZZ1int = ZZ1[,fv,]	
			# }else{
				# ZZ1int = apply(ZZ1[,fv,],c(1,3),sum)							
			# }
			# ZZ1 = array(ZZ1[,-fv,],dimz)
			# if(l==2) ZZ1int = matrix(ZZ1int,1,length(ZZ1int))				
			# sc_ga = est_multi_glob(YY,ZZ1,ltype,be=ga,Int=ZZ1int,only_sc=TRUE)$sc
		# }
		# ZZ = ZZ0
    	# for(j in 1:J){
			# ind = (refitem==j)
	    	# ZZ[,,ind] = ZZ[,,ind]*gac[j]
		# }
# #####		
#	}
	sc_par = est_multi_glob_gen(YY,ZZ,ltype,be=par,only_sc=TRUE)$sc
# Update piv
	if(k1==1) sc_de1 = NULL
	else sc_de1 = est_multi_glob(V%*%t(Aggr1),XX1dis,logit_cov,Xlabel,de1,only_sc=TRUE)$sc
	if(k2==1) sc_de2 = NULL
	else sc_de2 = est_multi_glob(V%*%t(Aggr2),XX2dis,logit_cov,Xlabel,de2,only_sc=TRUE)$sc
# output
	sc = c(sc_de1,sc_de2,sc_par,sc_ga1,sc_ga2)
	out = list(lk=lk,sc=sc)
	
}