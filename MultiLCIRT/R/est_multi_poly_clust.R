est_multi_poly_clust <- function(S,kU,kV,W=NULL,X=NULL,clust,start=0,link=0,disc=0,difl=0,
                           multi=1:J,piv=NULL,Phi=NULL,gac=NULL,DeU=NULL,DeV=NULL,fort=FALSE,tol=10^-10,
                           disp=FALSE,output=FALSE){

# With kV=1
	if(kV==1){
	  cat("fit only independence model\n")
	  kU = 1; link = 0; disc = 0; difl = 0; multi = 1:dim(S)[2]
	  X = NULL
	}
	if(max(S,na.rm=TRUE)==1 & difl!=0){
	  cat("with binary data put difl=0\n")
	  difl = 0		
	}
# Preliminaries
# check problems with input
    cov1 = !is.null(W); if(kU==1) cov1 = FALSE
    cov2 = !is.null(X); if(kV==1) cov2 = FALSE
    if(cov1){
	    	W = as.matrix(W)
    		namesW = colnames(W)
    }else{
    		namesW = NULL
    }
    if(cov2){
	    	X = as.matrix(X)
    		namesX = colnames(X)
	}else{
		namesX = NULL
	}
    miss = any(is.na(S))
	ns = nrow(S); J = ncol(S)
	nclust = max(clust)
    if(miss){
    	cat("Missing data in the dataset, units and items without responses are removed\n")
    	ind = which(apply(is.na(S),1,all))
    	if(length(ind)>0){
        	S = S[-ind,]
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
    if(cov1){
    	# WWdis has dimension equal to distinct W x kU
	    	ncov1 = ncol(W)
    		if(is.null(ncov1)){W = as.matrix(W); ncov1 = 1}
    		out = aggr_data(W,fort=fort)
    		Wdis = as.matrix(out$data_dis); Wlabel = out$label; Wndis = max(out$label)
    		WWdis = array(0,c(kU-1,(kU-1)*(1+ncov1),Wndis))
    		II = diag(kU-1)
   		for(i in 1:Wndis) WWdis[,,i] = II%x%t(c(1,Wdis[i,]))
   	}else{
   		ncov1 = 0
    		WWdis = array(diag(kU-1),c(kU-1,kU-1,1))
   		Wlabel = rep(1,nclust)
    }
    if(cov2){
# XXdis has dimension equal to distinct X x kV
	    	ncov2 = ncol(X)	    	
		if(is.null(ncov2)){X = as.matrix(X); ncov2 = 1}
    		out = aggr_data(X,fort=fort)
    		Xdis = as.matrix(out$data_dis); Xlabel = out$label; Xndis = max(out$label)
    		XXdis = array(0,c(kV-1,(kV-1)*(kU+ncov2),kU*Xndis))
    		II = diag(kV-1)
    		j = 0
    		for(u in 1:kU){
    			uv = rep(0,kU); uv[u] = 1
    			for(i in 1:Xndis){
    				j = j+1
    				XXdis[,,j] = II%x%t(c(uv,Xdis[i,]))
    			}
    			if(u>1) Xlabel = c(Xlabel,max(Xlabel)+out$label)
    		}
    }else{
	    	ncov2 = 0
	    	if(kV==1){
			XXdis = NULL
		}else{
			Xlabel = rep(1,ns)
    		XXdis = array(0,c(kV-1,(kV-1)*kU,kU))
    		II = diag(kV-1)
    		j = 0
    		for(u in 1:kU){
    			uv = rep(0,kU); uv[u] = 1
   				j = j+1
   				XXdis[,,j] = II%x%t(uv)
    			if(u>1) Xlabel = c(Xlabel,max(Xlabel)+rep(1,ns))
    		}			
		}
    }
# about models
	if(link==0){
		rm = J
		ltype = "m"
		ZZ = array(0,c(l-1,J*kV*(l-1),J*kV)); cont = 0
		for(c in 1:kV){
			u1 = matrix(0,1,kV); u1[c] = 1
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
#		DeU = matrix(0,J,rm)
		if(rm==1){
			DeU = 1
			fv = multi[1]
		}else{
			for(r in 1:rm){
				ind = multi[r,]
				ind = ind[ind>0]
#				DeU[ind,r] = 1      
			}
			fv = multi[,1]     # list of constrained items
		}
		fve = (fv-1)*(l-1)+1
		indga = 1:J; indga = indga[-fv]
		indth = 1:(kV*rm)
		if(difl==0){
			indbe = kV*rm+(1:(J*(l-1)-rm))
			indbec = 1:(J*(l-1)); indbec = indbec[-fve]
		}else{
			indbe = kV*rm+(1:(J-rm+rm*(l-2)));
			indbec = 1:J; indbec = indbec[-fv]
		}
# find non redundant response configurations for each item
		conf1 = NULL; conf2 = NULL
		for(j in 1:J) for(h in 0:l){
			if(any(S[,j]==h)){
				conf1 = c(conf1,j)
				conf2 = c(conf2,h)
			}
		}
		nconf = length(conf1)
# abililities for each item
		if(rm==1) abils=rep(1,J) else{
			abils = rep(0,J)
			for(h in 1:rm){
				ind = multi[h,]; ind = ind[ind>0]
				abils[ind] = h
			}
		}
# design matrix
		if(difl==0){
			ZZ = array(0,c(l-1,kV*rm+(l-1)*J-rm,J*kV))
		}else{
			if(difl==1) ZZ = array(0,c(l-1,kV*rm+J-rm+rm*(l-2),J*kV))
		}
		cont = 0; refitem = matrix(0,J*kV,1)       # reference item of that design matrix
		for(c in 1:kV){
			u1 = matrix(0,1,kV); u1[c] = 1
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
	    confe2 = rep(conf2,kV)   # response configuration
		confe1 = conf1          # corresponding design matrix in ZZ
		for(c in 2:kV) confe1 = c(confe1,max(confe1)+conf1)
	}
	ZZ0 = ZZ
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
		if(cov1){
			de1 = rep(0,(ncov1+1)*(kU-1))
		    la = NULL			
		}else{
			de1 = NULL; la = rep(1/kU,kU)
		}		
		if(cov2 | kV>1){
			de2 = rep(c(1:kU,rep(0,ncov2)),kV-1)
		    piv = NULL
	    }else{
			de2 = NULL; piv = rep(1,kV)/kV
		}
		if(kV==1) grid = 0 else grid = seq(-kV,kV,2*kV/(kV-1))
		Phi = array(0,c(l,J,kV)) # probability every response
		for(j in 1:J){
			dist = rep(0,l)
			for(y in 0:(l-1)) dist[y+1] = (sum(S[,j]==y)+0.5)/ns
			eta = Co%*%log(Ma%*%dist)
			for(c in 1:kV) Phi[,j,c] = inv_glob(eta+grid[c])$p
		}
	}
	if(start == 1){
		if(cov1){
			de1 = rnorm((kU-1)*(ncov1+1))/rep(c(1,apply(W,2,sd)),kU-1)
			la = NULL
		}else{
			de1 = NULL; la = runif(kU); la = la/sum(la)
		}
		if(cov2 || kV>1){
			de2 = rnorm((kV-1)*(ncov2+kU))/rep(c(rep(1,kU),apply(X,2,sd)),(kV-1))
			piv = NULL
		}else{
			de2 = NULL; piv = runif(kV); piv = piv/sum(piv)
		}
		Phi = array(runif(l*J*kV),c(l,J,kV))
		for(c in 1:kV) for(j in 1:J) Phi[,j,c] = Phi[,j,c]/sum(Phi[,j,c]);
	}
	if(start==2){
		de1 = as.vector(DeU)  
		de2 = as.vector(DeV)
	}  
	if(link==0){
		ga = NULL
	}else{
		if (start==0 || start==1) ga = rep(1,J-rm)
		else ga = gac[-fv]
	}
# Compute log-likelihood
	Psi = matrix(1,ns,kV) # probability observed response
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
	if(cov1){
	    out = prob_multi_glob(WWdis,"m",de1,Wlabel)
	    Ladis = out$Pdis; La = out$P
	}else{
		La = rep(1,nclust)%o%la
	}
	if(cov2 | kV>1){
	    out = prob_multi_glob(XXdis,"m",de2,Xlabel)
	    Pdis = out$Pdis; Piv = out$P
	 	Piv = array(t(Piv),c(kV,ns,kU))
	 	Piv = aperm(Piv,c(2,1,3))
	}else{
		Piv = rep(1,ns)%o%piv
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
	cat("*-------------------------------------------------------------------------------*\n")
	if(link==1 || link==2){
		cat(c("Multilevel model with multidimensional structure\n"))
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
    				cat("------------|-------------|-------------|-------------|-------------|\n")
	    			cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |\n");
    				cat("------------|-------------|-------------|-------------|-------------|\n")			
		}else{
			if(disc==0 || length(ga)==0){
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
	    			cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n");
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			}else if(disc==1){
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
				cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n");
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");

			}
			cat(sprintf("%11g",c(0,kU*10+kV,link,lk)),"\n",sep=" | ")
		}
	}
 	it = 0; lko = lk-10^10; dis = 0; par = NULL; dga = NULL; lkv = NULL
# Iterate until convergence
	while(((abs(lk-lko)/abs(lko)>tol) && it<10^4) || it<2){
		it = it+1
		paro = par; gao = ga; pivo = piv; de2o = de2; lko = lk
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
# Update item and ability parameters
		YY = matrix(0,J*kV,l)
		count = 0
		for(c in 1:kV) for(j in 1:J){
			count = count+1
			for(y in 1:l){
				ind = (S[,j]==(y-1))
				if(miss) YY[count,y] = sum(V[ind,c]*R[ind,j]) else YY[count,y] = sum(V[ind,c])			
			}
		}
		if(link==0){  # LC model
			cont = 0
			for(c in 1:kV) for(j in 1:J){
				cont = cont+1
				Phi[,j,c] = YY[cont,]/sum(YY[cont,])
			}
		}else{          
# other models
			if(disc==1){
				if(it>1 & rm<J){
					ZZ1 = array(0,c(l-1,J,J*kV))
					count = 0
					for(c in 1:kV) for(j in 1:J){
						count = count+1;
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
					ga = est_multi_glob(YY,ZZ1,ltype,be=ga,Int=ZZ1int)$be
				}
				gac = rep(1,J); gac[indga] = ga
				ZZ = ZZ0
    				for(j in 1:J){
	    				ind = (refitem==j)
		    			ZZ[,,ind] = ZZ[,,ind]*gac[j]
			    }
			}
			out = est_multi_glob(YY,ZZ,ltype,be=par)						
	    	par = out$be; P = out$P		    		
			Phi = array(t(P),c(l,J,kV))
    	}
# Update la
		if(cov1){
  		    out = est_multi_glob(Vclust,WWdis,"m",Wlabel,de1)
			de1 = out$be; Ladis = out$Pdis; La = out$P
		}else{
			la = colSums(Vclust)/nclust
			La = rep(1,nclust)%o%la
		}
# Update piv
        if(cov2 | kV>1){
  		    out = est_multi_glob(Vcomp,XXdis,"m",Xlabel,de2)
			de2 = out$be; Pdis = out$Pdis; Piv = out$P
		 	Piv = array(t(Piv),c(kV,ns,kU))
		 	Piv = aperm(Piv,c(2,1,3))
        }else{
            piv = sV/ns
            Piv = rep(1,ns)%o%piv
       }
# Compute log-likelihood
		Psi = matrix(1,ns,kV)
		if(miss){
			for(j in 1:J) for(c in 1:kV)	Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
		}else{
			if(fort){
                o = .Fortran("lk_obs",J,as.integer(kV),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
                Psi = o$Psi
            }else{
    			for(j in 1:J) for(c in 1:kV)	Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
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
		if(it>1 & link>0) dis = max(c(abs(par-paro),abs(ga-gao),abs(piv-pivo)))
		if(disp){
			if(it/10==floor(it/10)){
				if(link==0){
					cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko)),"\n",sep=" | ")				
				}else{
					if(disc==0 || length(ga)==0) cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ") else{
					if(disc==1) cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
					}
				}			
#t2 = proc.time()-t0
#print(c(t1[1],t2[1]))
			}
		}
		lkv = c(lkv,lk)
	}
	if(disp){
		if(it/10>floor(it/10)){
			if(link==0){
				cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko)),"\n",sep=" | ")				
			}else{
				if(disc==0 || length(ga)==0){
					cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
				}else if(disc==1){
					cat(sprintf("%11g",c(it,kU*10+kV,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
				}
			}
		}
		if(link==0){
	   		cat("------------|-------------|-------------|-------------|-------------|\n")
		}else{
			if(disc==0 || length(ga)==0){
		   		cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			}else if(disc==1){
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			}
		}
	}
# Compute number of parameters  
	if(link == 0){
		np = kV*J*(l-1)
    }else if(link==1 || link==2){
    		np = kV*rm+disc*(J-rm)
		if(difl==0) np = np+(l-1)*J-rm
		else if(difl==1) np = np+J-rm+l-2
	}
	np = np+(kV-1)*(ncov2+kU)
	np = np+(kU-1)*(ncov1+1)
# extract parameter estimates  
	aic = -2*lk+2*np;
	bic = -2*lk+np*log(ns);
	if(link==0){
		Th = NULL; Bec = NULL; gac = NULL; fv = NULL
	}else if(link==1 || link==2){
		th = par[indth]; be = par[indbe]
    		if(difl==0){
      		bec = rep(0,J*(l-1)); bec[indbec] = be
      		Bec = t(matrix(bec,l-1,J))
    		}else{
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
		Th = matrix(th,rm,kV)
    		names2 = NULL
    		for(c in 1:kV){names2 = c(names2,paste("Class",c))}
    		rownames(Th) = names
    		colnames(Th) = names2
  	}
# posterior probabilities
	if(kV==1){
		Vclust = matrix(1,nclust,1)
		Vind = matrix(1,ns,1)
	}else{
		Tmp = exp(lPj-mlPj)
		Vclust = Tmp*(1/rowSums(Tmp))
		Vcomp = NULL; Vind = 0
		Tmp1 = NULL
		for(u in 1:kU){
			Tmp = Psi*Piv[,,u]
			Tmp = Tmp*(Vclust[clust,u]/rowSums(Tmp))
			Vind = Vind+Tmp
		}
	}
# regression coefficients
	if(!cov1 & kU>1) de1 = log(la[2:kU]/la[1])
	if(kU>1){
		DeU = matrix(de1,ncov1+1,kU-1)
	}else{
		de1 = DeU = NULL
	}
	if(is.null(namesW) & ncov1>0) namesW = c("intercept",paste("W.",1:ncov1,sep=""))	
	else namesW = c("intercept",namesW)
	if(kU>1) dimnames(DeU) = list(namesW,logit=2:kU)
	if(cov2 | kV>1){
		DeV = matrix(de2,ncov2+kU,kV-1)
		if(is.null(namesX) & ncov2>0){
      		namesX = c(paste("intercept.",1:kU,sep=""),paste("X.",1:ncov2,sep=""))	
		}else{
			namesX = c(paste("intercept.",1:kU,sep=""),namesX)
	  	}
	    	dimnames(DeV) = list(namesX,logit=2:kV)
	}else{
  		DeV = NULL
  	}
	dimnames(Phi) = list(category=0:(l-1),item=1:J,class=1:kV)
	if(link==0){
  		par = NULL
  		for(c in 1:kV) for(j in 1:J) par = c(par,log(Phi[-1,j,c]/Phi[1,j,c]))
	}	
	if(output){
		lde1 = length(de1); lde2 = length(de2); lpar = length(par); lga = 0
		par_comp = c(de1,de2,par)
		if(disc==1){
			lga = length(ga)
			par_comp = c(par_comp,ga)
		}
		if(disp){
			cat("computation of derivatives\n")
			cat(length(par_comp),"parameters\n")
		}
  		out = lk_obs_score_clust(par_comp,lde1,lde2,lpar,lga,S,R,kU,kV,rm,l,J,fv,link,disc,indga,refitem,
                           miss,ltype,WWdis,Wlabel,XXdis,Xlabel,ZZ0,clust,fort)
		scn = rep(0,length(par_comp)); Jn = NULL
		for(j in 1:length(par_comp)){
			par_comp1 = par_comp; par_comp1[j] = par_comp1[j]+10^-6
			out1 = lk_obs_score_clust(par_comp1,lde1,lde2,lpar,lga,S,R,kU,kV,rm,l,J,fv,link,disc,indga,refitem,
                           miss,ltype,WWdis,Wlabel,XXdis,Xlabel,ZZ0,clust,fort)
			scn[j] = (out1$lk-lk)*10^6	
			Jn = cbind(Jn,(out1$sc-out$sc)*10^6)
			if(disp) if(j/10>floor(j/10)) cat(".") else cat(round(j/10))
	  		if(disp) if(j/100==floor(j/100)) cat("\n")
  		}
  		if(disp) cat("\n")
  		Jn = -(Jn+t(Jn))/2
  		Vn = ginv(Jn)
  		se = sqrt(abs(diag(Vn)))
  		sede1 = se[1:lde1]
  		sede2 = se[lde1+(1:lde2)]
  		separ = se[lde1+lde2+(1:lpar)]
  		if(disc==1) sega = se[lde1+lde2+lpar+(1:lga)] else sega = NULL
		if(kU>1){
			seDeU = matrix(sede1,ncov1+1,kU-1) 
			dimnames(seDeU) = list(namesW,logit=2:kU) 
		}else{
			seDeU = NULL
		}
		if(kV>1){
			seDeV = matrix(sede2,ncov2+kU,kV-1)
		    dimnames(seDeV) = list(namesX,logit=2:kV)
		}else{
			seDeV = NULL
		}
		# print(c(lk,out$lk,lk-out$lk))
  		# print(cbind(out$sc,scn,out$sc-scn))
	}  
  if(output){
    out = list(piv=piv,Th=Th,Bec=Bec,gac=gac,fv=fv,Phi=Phi,DeU=DeU,DeV=DeV,La=La,Piv=Piv,
           Vclust=Vclust,Vind=Vind,lk=lk,np=np,aic=aic,bic=bic,lkv=lkv,
           seDeU=seDeU,seDeV=seDeV,separ=separ,sega=sega,Vn=Vn,call=match.call())
  }else{
    out = list(piv=piv,Th=Th,Bec=Bec,gac=gac,fv=fv,DeU=DeU,DeV=DeV,La=La,Phi=Phi,
           lk=lk,np=np,aic=aic,bic=bic,call=match.call())
  }
  class(out) = "est_multi_poly_clust"
  out

}