est_multi_poly <- function(S,yv=rep(1,ns),k,X=NULL,start=0,link=0,disc=0,difl=0,
                           multi=NULL,piv=NULL,Phi=NULL,gac=NULL,De=NULL,fort=FALSE,tol=10^-10,
                           disp=FALSE,output=FALSE,out_se=FALSE,glob=FALSE){

#        [piv,Th,Bec,gac,fv,Phi,Pp,lk,np,aic,bic] = est_multi_poly(S,yv,k,start,link,disc,difl,multi,piv,Th,bec,gac)
#
# Fit Latent Class model and some restricted versions with k classes for ordinal (NA for missing data)
# 
# S    : matrix of available configurations
# X    : matrix of corresponding covariates affecting the ability
# yv   : frequencies of the available configurations
# k    : number of latent classes
# X    : matrix of covariates for the multinomial logit on the class weights
# start: type of startgine values (0 = deterministic, 1 = random, 2 = external input)
# link : type of link (0 = LC, 1 = GRM, 2 = PCM)
# disc : discriminating indices (0 = constrained, 1 = free)
# difl : difficulty levels (0 = free, 1 = additive decomposition)
# lk   : maximum log-likelihood
# piv  : weights of the latent classes
# Phi  : conditional distributions given the latent classes
# np   : number of free parameters
# bic  : Bayesian information criterion
# Th,be,ga : parameters for the model 4 (Th=Psi if start==3)
# fv   : list of items constrained
# fort : T for using fortran code for covariates, F using R code only
# tol  : relative tolerance level for convergence
# disp : TRUE for displying the likelihood evolution step by step
# output : to return additional outputs (Phi,Pp,Piv)
# out_se : to return standard errors
# glob : for global parametrization of the prior probabilities
#
# OUTPUT:
# ent  : entropy

# With k=1
	if(k==1){
	  cat("fit only for LC model with no other input\n")
	  link = 0; disc = 0; difl = 0#; multi = 1:dim(S)[2]
	  output = FALSE
	  X = NULL
	}
	if(max(S,na.rm=TRUE)==1 & difl!=0){
	  cat("with binary data put difl=0\n")
	  difl = 0		
	}
# Preliminaries
# check problems with input
    cov = !is.null(X)
    if(cov){
    	X = as.matrix(X)
   		namesX = colnames(X)
   		if(glob) logit_cov = "g" else logit_cov = "m"
    }else{
    	logit_cov = "m"
    }
    miss = any(is.na(S))
	ns = nrow(S); J = ncol(S)
    if(miss){
    	cat("Missing data in the dataset, units and items without responses are removed\n")
    	ind = which(apply(is.na(S),1,all))
    	if(length(ind)>0){
        	S = S[-ind,]; yv = yv[-ind]
        	if(!is.null(X)) X = as.matrix(X[-ind,])
        	ind = which(apply(is.na(S),2,all))
        	if(length(ind)>0){
        	    S = S[,-ind]
	            miss = any(is.na(S))
	        }
	    }
    }
    if(miss){R=1*(!is.na(S)); S[is.na(S)]=0}
	l = max(S)+1
	ns = nrow(S); J = ncol(S)
	n = sum(yv)
	if(link==0){
		 if(!is.null(multi)) warning("input multi not used")
	}else{
		if(is.null(multi)) multi = 1:J
	}
# check number of response categories
	flag = FALSE
	for(j in 1:J) if(length(table(S[,j]))<l) flag = TRUE
	if(flag){
		cat("|------------------------------------------ WARNING -----------------------------------------|\n")
	    cat("| Response variables must have the same number of categories starting with 0                 |\n")
	    cat("| If for some category the frequency is 0, use only LC models (link=0) or RS models (difl=1) |\n")	
		cat("|--------------------------------------------------------------------------------------------|\n\n")
	}
# checks about the covariates
    if(cov){
    	ncov = ncol(X)
    	out = aggr_data(X,fort=fort)
    	Xdis = as.matrix(out$data_dis); Xlabel = out$label; Xndis = max(out$label)
    	if(glob){
	    	XXdis = array(0,c(k-1,k-1+ncov,Xndis))
    		for(i in 1:Xndis) XXdis[,,i] = cbind(diag(k-1),rep(1,k-1)%o%Xdis[i,])    		
    	}else{
	    	XXdis = array(0,c(k-1,(k-1)*(ncov+1),Xndis))
    		if(k==2) II = 1 else II = diag(k-1)
    		for(i in 1:Xndis) XXdis[,,i] = II%x%t(c(1,Xdis[i,]))
    	}
    }else{
    	ncov = 0
    	XXdis = array(diag(k-1),c(k-1,k-1,1)); Xlabel = rep(1,ns)
    }
# about models
	if(link==0){
		rm = J
		if(glob) ltype = "g" else ltype = "m"
		ZZ = array(0,c(l-1,J*k*(l-1),J*k)); cont = 0
		for(c in 1:k){
			u1 = matrix(0,1,k); u1[c] = 1
			for(j in 1:J){
				cont = cont+1
				u2 = matrix(0,1,J); u2[j] = 1
				ZZ[,,cont] = (u1%x%u2)%x%diag(l-1)
			}
		}
	}
	if(link==1) ltype = "g" else if(link==2) ltype = "l"
	if(link == 1 || link == 2){
		if(is.vector(multi)) rm = 1
		else rm = nrow(multi)
		De1 = matrix(0,J,rm)
		if(rm==1){
			De1 = 1
			fv = multi[1]
		}else{
			for(r in 1:rm){
				ind = multi[r,]
				ind = ind[ind>0]
				De1[ind,r] = 1      
			}
			fv = multi[,1]     # list of constrained items
		}
		fve = (fv-1)*(l-1)+1
		indga = 1:J; indga = indga[-fv]
		indth = 1:(k*rm)
		if(difl==0){
			indbe = k*rm+(1:(J*(l-1)-rm))
			indbec = 1:(J*(l-1)); indbec = indbec[-fve]
		}else{
			indbe = k*rm+(1:(J-rm+rm*(l-2)));
			indbec = 1:J; indbec = indbec[-fv]
		}
# abililities for each item
		if(rm==1) abils=rep(1,J) else{
			abils = rep(0,J)
			for(h in 1:rm){
				ind = multi[h,]; ind = ind[ind>0]
				abils[ind] = h
			}
		}
# design matrix
		if(difl==0) ZZ = array(0,c(l-1,k*rm+(l-1)*J-rm,J*k))
		if(difl==1) ZZ = array(0,c(l-1,k*rm+J-rm+rm*(l-2),J*k))
		cont = 0; refitem = matrix(0,J*k,1)       # reference item of that design matrix
		for(c in 1:k){
			u1 = matrix(0,1,k); u1[c] = 1
			for(j in 1:J){
				u2 = matrix(0,1,rm); u2[abils[j]] = 1
				v = matrix(0,1,J); v[j] = 1
				cont = cont+1
				if(difl==0){
					Te = v%x%diag(l-1)
					Te = matrix(Te[,-fve],l-1,dim(Te)[2]-length(fve))
				}else if(difl==1){
					IS = diag(l-1); IS = IS[,-1]
					Te = cbind(v%x%rep(1,l-1),u2%x%IS)
					Te = Te[,-fv]
				}				
				ZZ[,,cont] = cbind(rep(1,l-1)%*%(u1%x%u2),-Te)
				refitem[cont] = j
			}
		}
	}
	ZZ0 = ZZ
# inequalities for ability
	if(glob){
		II = diag(k-1); II = cbind(0,II)-cbind(II,0)
		if(rm>1) II = II%x%matrix(c(1,rep(0,rm-1)),1,rm)
		Dis = cbind(II,matrix(0,k-1,dim(ZZ)[2]-k*rm))
	}
	else Dis = NULL
# When there is just 1 latent class
#   if k == 1,
#     piv = 1;
#     P = zeros(2,J);
#     for j in 1:J,
#       for jb = 0:1,
#         ind = which(S[,j]==jb);
#         P(jb+1,j) =  sum(yv(ind))/n;
#       end
#     end
#     Psi = P;
#     psi = ones(ns,1);
#     for j in 1:J,
#       psi = psi.*P(S[,j]+1,j);
#     end
#     lk = yv"*log(psi);
#     np = J;
#     aic = -2*lk+2*np;
#     bic = -2*lk+np*log(n);
#     bec = NULL; gac = NULL;
#     Pp = ones(ns,1);
#     Th = NULL;
#     return
#   end
	out = matr_glob(l); Co = out$Co; Ma = out$Ma
# Starting values
	if(start == 0){
		if(cov){
		    de = NULL; Piv = matrix(1/k,ns,k); piv = NULL				
	    }else{
			be = NULL; de = NULL; piv = rep(1,k)/k
		} # latent class probabilities
		if(k==1) grid = 0 else grid = seq(-k,k,2*k/(k-1))
		Phi = array(0,c(l,J,k)) # probability every response
		for(j in 1:J){
			dist = rep(0,l)
			for(y in 0:(l-1)) dist[y+1] = (sum(yv[S[,j]==y])+0.5)/n
			eta = Co%*%log(Ma%*%dist)
			for(c in 1:k) Phi[,j,c] = inv_glob(eta+grid[c])$p
		}
	}
	if(start == 1){
		if(cov){
			if(glob){
			    de = NULL; Piv = matrix(runif(ns*k),ns,k)
			    Piv = Piv*(1/rowSums(Piv))
			    piv = NULL
			}else{
				de = rnorm((k-1)*(ncov+1))/rep(c(1,apply(X,2,sd)),(k-1))
				Piv = prob_multi_glob(XXdis,logit_cov,de,Xlabel)$P
				piv = NULL
			}
		}else{
			piv = runif(k)
			piv = piv/sum(piv)
		}
		Phi = array(runif(l*J*k),c(l,J,k))
		for(c in 1:k) for(j in 1:J) Phi[,j,c] = Phi[,j,c]/sum(Phi[,j,c])
		if(glob){
			if(runif(1)>0.5) for(j in 1:J){
				mPhi = (0:(l-1))%*%Phi[,j,]
				ind = order(mPhi)
				Phi[,j,] = Phi[,j,ind]			
			}
		}
	}
	if(start==2) de = as.vector(De)  
	if(link==0){
		ga = NULL
	}else{
		if (start==0 || start==1) ga = rep(1,J-rm)
		else ga = gac[-fv]
	}
# Compute log-likelihood
	Psi = matrix(1,ns,k) # probability observed response
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
	if(cov){
		if(start==2) Piv = prob_multi_glob(XXdis,ltype,de,Xlabel)$P
	}else{
		Piv = rep(1,ns)%o%piv
	}
	if(k==1) Pj = Psi else Pj = Psi*Piv
	pm = rowSums(Pj)
	lk = sum(yv*log(pm))
	cat("*-------------------------------------------------------------------------------*\n")
	if(link==1 || link==2){
		cat(c("Model with multidimensional structure\n"))
		names = NULL
		for(j in 1:rm){names = c(names,paste("Dimension",j))}
		multi_out = as.matrix(multi)
		if(rm == 1) multi_out = t(multi_out)
		rownames(multi_out) = names
		print(multi_out)
	}
	cat(c("Link of type =                 ",link,"\n"))
	cat(c("Discrimination index =         ",disc,"\n"))
	cat(c("Constraints on the difficulty =",difl,"\n"))
	cat(c("Type of initialization =       ",start,"\n"))
	cat("*-------------------------------------------------------------------------------*\n")		
	if(disp){
		if(link==0){
    			cat("------------|-------------|-------------|-------------|-------------|\n");
    			cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |\n");
    			cat("------------|-------------|-------------|-------------|-------------|\n");
    		}else{
			if(disc==0 || length(ga)==0){
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    				cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n");
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
			}
			if(disc==1){
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
				cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n");
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
			}
		}
		cat(sprintf("%11g",c(0,k,link,lk)),"\n",sep=" | ")
	}
 	it = 0; lko = lk-10^10; dis = 0; par = NULL; dga = NULL; lkv = NULL
# Iterate until convergence
	while(((abs(lk-lko)/abs(lko)>tol) && it<10^4) || it<2){
#t0 = proc.time()
		it = it+1
		paro = par; gao = ga; pivo = piv; deo = de; lko = lk
# ---- E-step ----
		V = ((yv/pm)%o%rep(1,k))*Piv*Psi; sV = colSums(V)
#print(proc.time()-t0)
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
		if(link==0){  # LC model
			if(glob){
				out = est_multi_glob(YY,ZZ,ltype,be=par,Dis=Dis)						
		    	par = out$be; P = out$P		    		
				Phi = array(t(P),c(l,J,k))
			}else{
				cont = 0
				for(c in 1:k) for(j in 1:J){
					cont = cont+1
					Phi[,j,c] = YY[cont,]/sum(YY[cont,])
				}
			}
		}else{          # other models
#print(proc.time()-t0)			
			if(disc==1){
				if(it>1 & rm<J){
					ZZ1 = array(0,c(l-1,J,J*k))
					count = 0
					for(c in 1:k) for(j in 1:J){
						count = count+1;
						ZZ1[,j,count] = ZZ0[,,count]%*%par
					}
					dimz = dim(ZZ1)
					dimz[2] = dimz[2]-length(fv)
					if(rm==1){
						ZZ1int = ZZ1[,fv,]	
					}else{
						Tmp = array(ZZ1[,fv,],c(l-1,rm,J*k))
						ZZ1int = apply(Tmp,c(1,3),sum)							
					}
					ZZ1 = array(ZZ1[,-fv,],dimz)
					if(l==2) ZZ1int = matrix(ZZ1int,1,length(ZZ1int))
					ga = est_multi_glob(YY,ZZ1,ltype,be=ga,Int=ZZ1int)$be
				}
				gac = rep(1,J); gac[indga] = ga
				ZZ = ZZ0
    				for(j in 1:J){
	    				ind = (refitem==j)
		    			ZZ[,,ind] = ZZ[,,ind]*gac[j]
			    }
			}
			out = est_multi_glob(YY,ZZ,ltype,be=par,Dis=Dis)						
	    	par = out$be; P = out$P		    		
			Phi = array(t(P),c(l,J,k))
    	}
# Update piv
		if(cov){
			out = est_multi_glob(V,XXdis,logit_cov,Xlabel,de)
	    	de = out$be; Pdis = out$Pdis; Piv = out$P
		}else{
			piv = sV/n
			Piv = rep(1,ns)%o%piv
		} 
#print(proc.time()-t0)
# Compute log-likelihood
		Psi = matrix(1,ns,k);
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
		if(it>1 & link>0) dis = max(c(abs(par-paro),abs(ga-gao),abs(piv-pivo)))
		if(disp){
			if(it/10==floor(it/10)){
				if(link==0){
					cat(sprintf("%11g",c(it,k,link,lk,lk-lko)),"\n",sep=" | ")
				}else{
					if(disc==0 || length(ga)==0) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
					if(disc==1) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
				}
			}
		}
		lkv = c(lkv,lk)
	}
	if(disp){
		if(it/10>floor(it/10)){
			if(link==0){
				cat(sprintf("%11g",c(it,k,link,lk,lk-lko)),"\n",sep=" | ")
			}else{
				if(disc==0 || length(ga)==0) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
				if(disc==1) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
			}
		}
		if(link==0){
    			cat("------------|-------------|-------------|-------------|-------------|\n")
		}else{
			if(disc==0 || length(ga)==0) cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			if(disc==1) cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
		}
	}
# Compute number of parameters  
	if(link == 0){
	  np = k*J*(l-1)
      if(cov) np = np+(k-1)*(ncov+1) else np = np+k-1
    }else if(link==1 || link==2){
      np = k*rm+disc*(J-rm)
      if(cov){
      	if(glob){
	      	np = np+k-1+ncov      		
      	}else{
	      	np = np+(k-1)*(ncov+1)
	    }
      }else{
      	np = np+k-1
      }
		if(difl==0) np = np+(l-1)*J-rm
		else if(difl==1) np = np+J-rm+l-2
	}
# extract parameter estimates  
  aic = -2*lk+2*np
  bic = -2*lk+np*log(n)
  if(link==0){
    Th = NULL; Bec = NULL; gac = NULL; fv = NULL
  }
  else if(link==1 || link==2){
    th = par[indth]; be = par[indbe]
    if(difl==0){
      bec = rep(0,J*(l-1)); bec[indbec] = be
      Bec = t(matrix(bec,l-1,J))
    }
    else{
      bec1 = rep(0,J); bec1[indbec] = be[1:(J-rm)]
      if(rm==1){
	      bec2 = rep(0,l-1); bec2[2:(l-1)] = be[J-rm+(1:(l-2))]
      }else{
	      bec2 = matrix(0,l-1,rm); bec2[2:(l-1),] = be[J-rm+(1:(rm*(l-2)))]
    	  dimnames(bec2) = list(level=1:(l-1),dim=1:rm)
   	  }
      Bec = list(difficulties=bec1,cutoffs=bec2)
    }
    gac = rep(1,J); gac[indga] = ga
    Th = matrix(th,rm,k)
    names2 = NULL
    for(c in 1:k){names2 = c(names2,paste("Class",c))}
    rownames(Th) = names
    colnames(Th) = names2
  }
  Pp = ((1./pm)%o%rep(1,k))*Piv*Psi
  ent = -sum(V*log(pmax(Pp,10^-100)))
  if(cov){
		if(glob){
			De = matrix(de,ncov+k-1,1)
			names_cutoff = paste("cutoff",1:(k-1),sep="")
   	   		if(is.null(namesX)){
      			namesX = c(names_cutoff,paste("X",1:ncov,sep=""))	
      		}else{
				namesX = c(names_cutoff,namesX)
	  		}
	  		rownames(De) = namesX
		}else{
	 		De = matrix(de,ncov+1,k-1)
   	   		if(is.null(namesX)){
      			namesX = c("intercept",paste("X",1:ncov,sep=""))	
      		}else{
				namesX = c("intercept",namesX)
	  		}
		    dimnames(De) = list(namesX,logit=2:k)
	  	}
		piv = colMeans(Piv)
  }else De = NULL
	dimnames(Phi) = list(category=0:(l-1),item=1:J,class=1:k)
	if(link==0 & !glob){
		Phi = pmax(Phi,10^-50)
  		par = NULL
  		for(c in 1:k) for(j in 1:J) par = c(par,log(Phi[-1,j,c]/Phi[1,j,c]))
	}
	if(!cov){
		de = De = log(piv[-1]/piv[1])
	}
	if(out_se){
		lde = length(de); lpar = length(par); lga = 0
		par_comp = c(de,par)
		if(disc==1){
			lga = length(ga)
			par_comp = c(par_comp,ga)
		}
		if(disp){
			cat("computation of derivatives\n")
			cat(length(par_comp),"parameters\n")
		}		
  		out = lk_obs_score(par_comp,lde,lpar,lga,S,R,yv,k,rm,l,J,fv,link,disc,indga,glob,refitem,
                           miss,ltype,XXdis,Xlabel,ZZ0,fort)
		scn = rep(0,length(par_comp)); Jn = NULL
		for(j in 1:length(par_comp)){
			par_comp1 = par_comp; par_comp1[j] = par_comp1[j]+10^-6
			out1 = lk_obs_score(par_comp1,lde,lpar,lga,S,R,yv,k,rm,l,J,fv,link,disc,indga,glob,refitem,
   	                            miss,ltype,XXdis,Xlabel,ZZ0,fort)
			 scn[j] = (out1$lk-lk)*10^6	
			 Jn = cbind(Jn,(out1$sc-out$sc)*10^6)
			if(disp) if(j/10>floor(j/10)) cat(".") else cat(round(j/10))
	  		if(disp) if(j/100==floor(j/100)) cat("\n")		
  		}
  		if(disp) cat("\n")  		
  		Jn = (Jn+t(Jn))/2
  		Vn = ginv(Jn)
  		se = sqrt(abs(diag(Vn)))
  		sede = se[1:lde]
  		separ = se[lde+(1:lpar)]
  		if(disc==1) sega = se[lpar+lde+(1:lga)] else sega = NULL
		if(glob) seDe = matrix(sede,ncov+k-1,1)
		else seDe = matrix(sede,ncov+1,k-1)  	
		# print(c(lk,out$lk,lk-out$lk))
  		# print(cbind(out$sc,scn,out$sc-scn))
	}
    out = list(piv=piv,Th=Th,Bec=Bec,gac=gac,fv=fv,De=De,Phi=Phi,
           lk=lk,np=np,aic=aic,bic=bic,ent=ent,call=match.call())
    if(output){	out$Piv=Piv; out$Pp=Pp; out$lkv = lkv}
    if(out_se){out$seDe=seDe; out$separ=separ; out$sega=sega; out$Vn=Vn}
	class(out) = "est_multi_poly"
  	return(out)

}