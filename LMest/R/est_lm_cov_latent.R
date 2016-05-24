est_lm_cov_latent <-
function(S,X1,X2,yv=rep(1,nrow(S)),k,start=0,tol=10^-8,maxit=1000,param="multilogit",
                       Psi,Be,Ga,fort=TRUE,output=FALSE,out_se=FALSE,fixPsi=FALSE){

# Fit the LM model with individual covariates in the distribution of the latent process
#
# INPUT:
# S = array of available configurations (n x TT x r)
# X1 = matrix of covariates affecting the initial probabilities
# X2 = array of covariates affecting the transition probabilities
# Psi = conditional response probabilities
# Be = parameters on the initial probabilities (if start=2)
# Ga = parameters on the transition probabilities (if start=2)
# start = initialization (0 = deterministic, 1 = random, 2 = initial values in input)
# param = type of parametrization for the transition probabilities:
#         multilogit = standard multinomial logit for every row of the transition matrix
#         difflogit  = multinomial logit based on the difference between two sets of parameters
# fort   = fortran use (FALSE for not use fortran)
# output = to return additional output
# out_se  = TRUE for computing the information and standard errors
# fixPsi = TRUE if Psi is given in input and is not updated anymore 

# Preliminaries
    check_der = FALSE # to check score and info
    if(fort!=TRUE) fort = FALSE
   	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
	n = sum(yv)
  	if(length(sS)==2) r = 1
  	else r = sS[3]
  	if(min(S,na.rm=TRUE)>0){
  		cat("|------------------- WARNING -------------------|\n")
  		cat("|The first response category must be coded as 0 |\n")
  		cat("|-----------------------------------------------|\n")	
 	} 
  	if(is.data.frame(S)){
  		warning("Data frame not allowed for S")
  	}
  	if(ns!=length(yv)) stop("dimensions mismatch between S and yv")
   	if(any(is.na(X1)) | any(is.na(X2))) stop("missing data not allowed in X1 or X2")
  	miss = any(is.na(S))
	if(miss){
         cat("Missing data in the dataset, treated as missing at random\n")
         R = 1 * (!is.na(S))
         S[is.na(S)] = 0
	}else{
		R = NULL
	}
  	
	Sv = matrix(S,ns*TT,r)
	if(miss) Rv = matrix(R,ns*TT,r)
	
	if(r==1){
		if(is.matrix(S)) S = array(S,c(dim(S),1))
		b = max(S); mb = b; sb = b
	    Co = cbind(-diag(b),diag(b))
	 	Ma = cbind(lower.tri(matrix(1,b,b), diag = TRUE),rep(0,b))
  		Ma = rbind(Ma,1-Ma)
  	}else{
		b = rep(0,r)
		for(j in 1:r) b[j] = max(S[,,j])
		mb = max(b); sb = sum(b)
		Matr = vector("list",r)
		for(j in 1:r){
			Matr[[j]]$Co = cbind(-diag(b[j]),diag(b[j]))
		 	Ma = cbind(lower.tri(matrix(1,b[j],b[j]), diag = TRUE),rep(0,b[j]))
  			Matr[[j]]$Ma = rbind(Ma,1-Ma)
		}
	}
	th = NULL; sc = NULL
	J = NULL
  		
# Covariate structure and related matrices: initial probabilities
	if(is.vector(X1)) X1 = matrix(X1,ns,1)
	nc1 = dim(X1)[2] # number of covariates on the initial probabilities
	if(ns!= dim(X1)[1]) stop("dimension mismatch between S and X1")
 
	nameBe = colnames(X1)
	if(k == 2) GBe = as.matrix(c(0,1)) else{
		GBe = diag(k); GBe = GBe[,-1]
	}
	out = aggr_data(X1,fort=fort)
	Xdis = out$data_dis
	if(nc1==1) Xdis = matrix(Xdis,length(Xdis),1)
	Xlab = out$label
	Xndis = max(Xlab)
	XXdis = array(0,c(k,(k-1)*(nc1+1),Xndis))
	for(i in 1:Xndis){
		xdis = c(1,Xdis[i,])
		XXdis[,,i] = GBe%*%(diag(k-1)%x%t(xdis))
	}

# for the transition probabilities
	if(TT==2) X2 = array(X2,c(ns,1,dim(X2)[2]))
	if(is.matrix(X2)) X2 = array(X2,c(ns,TT-1,1))
    nc2 = dim(X2)[3] # number of covariates on the transition probabilities
    if(ns!= dim(X2)[1]) stop("dimension mismatch between S and X2")

    nameGa = colnames(aperm(X2,c(1,3,2)))
	Z = NULL
	for(t in 1:(TT-1)) Z = rbind(Z,X2[,t,])
	if(nc2==1) Z = as.vector(X2)
	out = aggr_data(Z,fort=fort); Zdis = out$data_dis; Zlab = out$label; Zndis = max(Zlab)
	if(nc2==1) Zdis=matrix(Zdis,length(Zdis),1)
	if(param=="multilogit"){
    	ZZdis = array(0,c(k,(k-1)*(nc2+1),Zndis,k))
	    for(h in 1:k){
		    if(k==2){
			    if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
		    }else{
			    GGa = diag(k); GGa = GGa[,-h]
		    }  		
		    for(i in 1:Zndis){
			    zdis = c(1,Zdis[i,])
			    ZZdis[,,i,h] = GGa%*%(diag(k-1)%x%t(zdis))
		    }
	    }
	 }else if(param=="difflogit"){
        Zlab = (((Zlab-1)*k)%x%rep(1,k))+rep(1,ns*(TT-1))%x%(1:k)
        ZZdis = array(0,c(k,k*(k-1)+(k-1)*nc2,Zndis*k))
        j = 0
		for(i in 1:Zndis){
            for(h in 1:k){
                j = j+1
                if(k==2){
			      if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
		        }else{
			        GGa = diag(k); GGa = GGa[,-h]
		        }  		
			    u = matrix(0,1,k); u[1,h] = 1
			    U = diag(k); U[,h] = U[,h]-1
			    U = U[,-1]
		        ZZdis[,,j] = cbind(u%x%GGa,U%x%t(Zdis[i,]))            
            }
	    }
    }
	
# for information matrix
  	if(out_se){
  	    Am = vector("list",r)
  	    for(j in 1:r) Am[[j]] = rbind(rep(0,b[j]),diag(b[j]))
  	}
# When there is just 1 latent class
  	if(k == 1){
		Piv = rep(1,n); Pi = 1
   	 	P = matrix(0,mb+1,r)
	    	for(t in 1:TT){
	      		for(j in 1:r){
		        	for(y in 0:b[j]){
		        		ind = which(S[,t,j]==y)
		  		        P[y+1,j] = P[y+1,j]+sum(yv[ind])
				}
	   		}
	    }
	    Psi = P/(n*TT)
	    pm = rep(1,ns)
	    if (miss) for(t in 1:TT) for(j in 1:r)  pm = pm*(Psi[S[,t,j]+1,j]*R[,t,j]+(1-R[,t,j]))
	    else for(t in 1:TT) for(j in 1:r)  pm = pm*Psi[S[,t,j]+1,j]
	    
	    lk = sum(yv*log(pm))
		if(r==1) np = k*mb*r else np = k*sum(b)
	    aic = -2*lk+np*2
	    bic = -2*lk+np*log(n)
		out = list(lk=lk,Piv=Piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=NULL,V=NULL,call=match.call())
		class(out)="LMlatent"
	    return(out)
  	}
    time = proc.time()
# Starting values: deterministic initialization
	if(start == 0){
		if(fixPsi==FALSE){
			P = matrix(NA,mb+1,r)
			for(j in 1:r) P[1:(b[j]+1),j] = 0
	        for(t in 1:TT) for(j in 1:r) for(y in 0:mb){
	        	ind = which(S[,t,j]==y)
		        if(miss) P[y+1,j] = P[y+1,j]+sum(t(R[ind,t,j])%*%yv[ind])
		        else P[y+1,j] = P[y+1,j]+sum(yv[ind])
	        }
	       	if(r==1) E = Co%*%log(Ma%*%P) else{
	       		E = matrix(NA,mb,r)
	    		for(j in 1:r){
	    			Co = Matr[[j]]$Co; Ma = Matr[[j]]$Ma
	    			E[1:b[j],j] = Co%*%log(Ma%*%P[1:(b[j]+1),j])
				}   		
	       	}
		  	Psi = array(NA,c(mb+1,k,r)); Eta = array(NA,c(mb,k,r))
	        grid = seq(-k,k,2*k/(k-1))
	        for(c in 1:k){
	         	for(j in 1:r){
		     		etac = E[1:b[j],j]+grid[c]
	          		Eta[1:b[j],c,j] = etac
	          		Psi[1:(b[j]+1),c,j] = invglob(etac)
	          	}
	       	 }  
       	 } 
# parameters on initial probabilities
       	be = array(0,(nc1+1)*(k-1))
       	out = prob_multilogit(XXdis,be,Xlab,fort)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
            Ga = matrix(0,(nc2+1)*(k-1),k)
            Ga[1+(0:(k-2))*(nc2+1),] = -log(10)
			PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,ns,TT))
			for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,ns,TT-1))
		    }
		}else if(param=="difflogit"){
             Ga = matrix(0,k*(k-1)+(k-1)*nc2)
             Ga[1:((h-1)*k)] = -log(10)
             PI = array(0,c(k,k,ns,TT))
             out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		     PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,ns,TT-1))
   		     PI = aperm(PI,c(2,1,3,4))
		}
  	}
# random initialization
  	if(start==1){
  		if(fixPsi==FALSE){
		    Psi = array(NA,c(mb+1,k,r)) 
		   # for(j in 1:r) Psi[1:(b[j]+1),,j] = 0
		    for(j in 1:r){
			    Psi[1:(b[j]+1),,j] = matrix(runif((b[j]+1)*k),b[j]+1,k) 
			    for(c in 1:k) Psi[1:(b[j]+1),c,j] = Psi[1:(b[j]+1),c,j]/sum(Psi[1:(b[j]+1),c,j])
		    }
	    }
# parameters on initial probabilities
       	be = c(rnorm(1),rep(0,nc1))
       	if(k>2) for(h in 2:(k-1)) be = c(be,rnorm(1),rep(0,nc1))
       	out = prob_multilogit(XXdis,be,Xlab,fort)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
    	#	Ga = matrix(-abs(rnorm((nc2+1)*(k-1),k)),(nc2+1)*(k-1),k)/2
    		Ga = matrix(0,(nc2+1)*(k-1),k)
    		Ga[1+(0:(k-2))*(nc2+1),] = -abs(rnorm((k-1)))   
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,ns,TT))
		    for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,ns,TT-1))
		   }
    	}else if(param=="difflogit"){
            Ga = c(-abs(rnorm(k*(k-1))),rep(0,(k-1)*nc2))
            PI = array(0,c(k,k,ns,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		    PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,ns,TT-1))
   		    PI = aperm(PI,c(2,1,3,4))
		}
	}
# initialization as input
    if(start==2){
# parameters on initial probabilities
       	be = as.vector(Be)
       	out = prob_multilogit(XXdis,be,Xlab,fort)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
        	if(is.list(Ga)) stop("invalid mode (list) for Ga")
    		Ga = matrix(Ga,(nc2+1)*(k-1),k)
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,ns,TT))
		    for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,ns,TT-1))
		    }
		}else if(param=="difflogit"){
			if(is.list(Ga)) Ga = c(as.vector(t(Ga[[1]])),as.vector(Ga[[2]]))
            if(length(Ga)!=k*(k-1)+(k-1)*nc2) stop("invalid dimensions for Ga")
            PI = array(0,c(k,k,ns,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		    PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,ns,TT-1))
   		    PI = aperm(PI,c(2,1,3,4))
		}
    }
  
###### standard EM #####
   	out = lk_comp_latent(S,R,yv,Piv,PI,Psi,k,fort=fort)    	
   	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  # 	if(is.nan(lk)) browser()
	it = 0; lko = lk-10^10; lkv = NULL
	par = c(as.vector(Piv),as.vector(PI),as.vector(Psi))
	if(any(is.na(par))) par = par[-which(is.na(par))]
	paro = par
# Iterate until convergence
# display output	
cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("      k     |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(k,start,0,lk)),"\n",sep=" | ")
   	#cat("",sprintf("%11g",c(0,lk)),"\n",sep=" | ")
	while((lk-lko)/abs(lk)>tol & it<maxit){
		Psi0 = Psi; Piv0 = Piv; PI0 = PI
		it = it+1
# ---- E-step ----
# Compute V and U
       	out = prob_post_cov(S,yv,Psi,Piv,PI,Phi,L,pv,fort=fort)
       	U = out$U; V = out$V      		
# If required store parameters
# ---- M-step ----
# Update Psi
       	if(fixPsi==FALSE){
    		Y1 = array(NA,c(mb+1,k,r))
    		for(j in 1:r) Y1[1:(b[j]+1)] = 0
	    	Vv = matrix(aperm(V,c(1,3,2)),ns*TT,k)
		    for(j in 1:r) for(jb in 0:b[j]) {
	    		ind = which(Sv[,j]==jb)
	    		if(length(ind)==1){
	    			if(miss) Y1[jb+1,,j] = Vv[ind,]*Rv[ind,j]
	    				else Y1[jb+1,,j] = Vv[ind,]
	    		}
	    		if(length(ind)>1){
	    				if(miss) Y1[jb+1,,j] = colSums(Vv[ind,]*Rv[ind,j])
	    				else Y1[jb+1,,j] = colSums(Vv[ind,])
	    		}

	    						
	    	}
	    	for(j in 1:r) for(c in 1:k){
				tmp = Y1[1:(b[j]+1),c,j]
				tmp = pmax(tmp/sum(tmp),10^-10)
				Psi[1:(b[j]+1),c,j] = tmp/sum(tmp)
			}
		}
# Update piv
		out = est_multilogit(V[,,1],XXdis,Xlab,be,Pivdis,fort=fort)
		be = out$be; Pivdis = out$Pdi; Piv = out$P
# Update Pi
		if(param=="multilogit"){
	    	for(h in 1:k){
		    	UU = NULL
		    	for(t in 2:TT) UU = rbind(UU,t(U[h,,,t]))
		    	out = est_multilogit(UU,ZZdis[,,,h],Zlab,Ga[,h],PIdis[,,h],fort=fort)
		    	PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,ns,TT-1)); Ga[,h] = out$be
	   	 	}
		}else if(param=="difflogit"){
		    Tmp = aperm(U[,,,2:TT],c(1,3,4,2))
		    Tmp = matrix(Tmp,ns*k*(TT-1),k)
           	out = est_multilogit(Tmp,ZZdis,Zlab,Ga,PIdis,fort=fort)
	    	PIdis = out$Pdis; Ga = out$be
	    	Tmp = array(out$P,c(k,ns,TT-1,k))
	    	PI[,,,2:TT] = aperm(Tmp,c(1,4,2,3)) 
        }
# Compute log-likelihood
   		paro = par; par = c(as.vector(Piv),as.vector(PI),as.vector(Psi))
   		if(any(is.na(par))) par = par[-which(is.na(par))]
   		lko = lk;
   		out = lk_comp_latent(S,R,yv,Piv,PI,Psi,k,fort=fort)
   		lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
# Display output
       	if(it/10 == floor(it/10)){
       		#cat("",sprintf("%11g",c(it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
       		cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
       	}
   		lkv = c(lkv,lk)
	}
if(it/10 > floor(it/10))  cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
#### compute infomation matrix ####
	if(out_se){
		dlPsi = array(NA,c(mb+1,k,r,k*sb))
		for(j in 1:r) dlPsi[1:(b[j]+1),,j,] = 0
   		count = 0
		for(c in 1:k) for(j in 1:r){
	  		ind = count+(1:b[j])
			temp = pmax(Psi[1:(b[j]+1),c,j],10^-50)
			dlPsi[1:(b[j]+1),c,j,ind] = (diag(b[j]+1)-rep(1,b[j]+1)%o%temp)%*%Am[[j]]
			count = count+b[j]
			th = c(th,log(temp[-1]/temp[1]))
		}
		dlPiv = array(0,c(ns,k,(1+nc1)*(k-1)))
		for(j in 1:Xndis){
	        temp = pmax(Pivdis[j,],10^-50)
	        Temp = (diag(k)-rep(1,k)%o%temp)%*%XXdis[,,j]
	    	for(i in which(Xlab==j)) dlPiv[i,,] = Temp
        }
        th = c(th,be)
       
        count = 0
        if(param=="multilogit"){
        	dlPI = array(0,c(k,k,ns*TT,(1+nc2)*(k-1)*k))
        	temp0 = rep(1,k); Temp0 = diag(k)
		    for(h in 1:k){
		    	ind = count+(1:((1+nc2)*(k-1)))
		    	for(j in 1:Zndis){
		    		temp = pmax(PIdis[j,,h],10^-50)
           	        Temp = (Temp0-temp0%o%temp)%*%ZZdis[,,j,h]
           	        for(i in which(Zlab==j)) dlPI[h,,ns+i,ind] = Temp
           	   	}
		    	count = count+((1+nc2)*(k-1))
                th = c(th,Ga[,h])		    
		    }
            dlPI = array(dlPI,c(k,k,ns,TT,(1+nc2)*(k-1)*k))
		}else if(param=="difflogit"){
		   dlPI = array(0,c(k,k*ns*TT,(k+nc2)*(k-1)))
		   temp0 = rep(1,k); Temp0 = diag(k)		   
		   for(j in 1:(Zndis*k)){	
		   		 temp = pmax(PIdis[j,],10^-50)
           	     Temp = (Temp0-temp0%o%temp)%*%ZZdis[,,j]
           	     for(i in which(Zlab==j)) dlPI[,k*ns+i,] = Temp         	 
		    	}
#		    Ga = c(as.vector(t(Ga[[1]])),as.vector(Ga[[2]]))
		    th = c(th,Ga)		    
            dlPI = array(dlPI,c(k,k,ns,TT,(k+nc2)*(k-1)))
            dlPI = aperm(dlPI,c(2,1,3,4,5))
		}
		
# Compute log-likelihood
		lk2 = lk
    	out = lk_comp_latent(S,R,yv,Piv,PI,Psi,k,der=TRUE,fort=fort,dlPsi=dlPsi,dlPiv=dlPiv,dlPI=dlPI)
    	sc = out$dlk; dlL = out$dlL; dlPhi = out$dlPhi; dlL2 = out$dlL2; dlpv = out$dlpv
    	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    	sc2 = sc; 
  #      if(is.nan(lk)) browser()
  		it = 0; lko = lk-10^10; lkv = NULL; dev = NULL
# backward recursion
        out = prob_post_cov(S,yv,Psi,Piv,PI,Phi,L,pv,der=TRUE,fort=fort,dlPhi=dlPhi,dlPiv,
   	                        dlPI=dlPI,dlL=dlL,dlL2=dlL2,dlpv=dlpv)
       	U = out$U; V = out$V; dlU = out$dlU; dlV = out$dlV
# ---- M-step ----
# score and info Psi
       	sc = NULL
   		Y1 = array(NA,c(mb+1,k,r))
   		for(j in 1:r) Y1[1:(b[j]+1),,j] = 0
    	Vv = matrix(aperm(V,c(1,3,2)),ns*TT,k)
	    for(j in 1:r) for(jb in 0:b[j]) {
	    	ind = which(Sv[,j]==jb)
	    	if(length(ind)==1){
	    		if(miss) Y1[jb+1,,j] = Vv[ind,]*Rv[ind,j]
	    			else Y1[jb+1,,j] = Vv[ind,]
	    	}
	    	if(length(ind)>1){
	    		if(miss) Y1[jb+1,,j] = colSums(Vv[ind,]*Rv[ind,j])
	    			else Y1[jb+1,,j] = colSums(Vv[ind,])
	    	}
	    					
	    }
	    for(c in 1:k) for(j in 1:r){
	    	sc = c(sc,t(Am[[j]])%*%(Y1[1:(b[j]+1),c,j]-sum(Y1[1:(b[j]+1),c,j])*Psi[1:(b[j]+1),c,j]))
	    	#Psi[,c,j] = Y1[,c,j]/sum(Y1[,c,j])
	    	tmp = Y1[1:(b[j]+1),c,j]
			tmp = pmax(tmp/sum(tmp),10^-10)
			Psi[1:(b[j]+1),c,j] = tmp/sum(tmp)
			
	        temp = pmax(Psi[1:(b[j]+1),c,j],10^-50)
            Op = diag(temp)-temp%o%temp
	        Temp  = sum(Y1[1:(b[j]+1),c,j])*t(Am[[j]])%*%Op%*%Am[[j]]
            if(j==1 & c==1) Fi = Temp else Fi = blkdiag(Fi,Temp)
	    }
# score and info piv
		out = est_multilogit(V[,,1],XXdis,Xlab,be,Pivdis,fort=fort,ex=TRUE)
		sc = c(sc,out$sc); Fi = blkdiag(Fi,out$Fi)
# score and info Pi
        if(param=="multilogit"){
		    for(h in 1:k){
			    UU = NULL
			    for(t in 2:TT) UU = rbind(UU,t(U[h,,,t]))
		    	out = est_multilogit(UU,ZZdis[,,,h],Zlab,Ga[,h],PIdis[,,h],fort=fort,ex=TRUE)
			    sc = c(sc,out$sc); Fi = blkdiag(Fi,out$Fi) 
	    	}
		}else if(param=="difflogit"){
	    	Tmp = aperm(U[,,,2:TT],c(1,3,4,2))
	    	Tmp = matrix(Tmp,ns*k*(TT-1),k)
           	out = est_multilogit(Tmp,ZZdis,Zlab,Ga,PIdis,fort=fort,ex=TRUE)
        	sc = c(sc,out$sc); Fi = blkdiag(Fi,out$Fi)
       	}
       	Fi = as.matrix(Fi)
# compute correction matrix for the information
   		nal = dim(dlPhi)[4]; nbe = dim(dlPiv)[3]; nga = dim(dlPI)[5]
		npar = nal+nbe+nga
		Cor = matrix(0,npar,npar)
		dY1 = array(NA,c(mb+1,k,r,npar))
		for(j in 1:r) dY1[1:(b[j]+1),,j,] = 0
		dV = array(V,c(ns,k,TT,npar))*dlV
    	dVv = array(aperm(dV,c(1,3,2,4)),c(ns*TT,k,nal+nbe+nga))
		for(j in 1:r) for(jb in 0:b[j]) {
			ind = which(Sv[,j]==jb)
			if(length(ind)==1){
				if(miss) for(h in 1:(nal+nbe+nga)) dY1[jb+1,,j,h] = dVv[ind,,h]*Rv[ind,j]
				else for(h in 1:(nal+nbe+nga)) dY1[jb+1,,j,h] = dVv[ind,,h]
			}
			if(length(ind)>1){
				if(miss) for(h in 1:(nal+nbe+nga)) dY1[jb+1,,j,h] = colSums(dVv[ind,,h]*Rv[ind,j])
				else for(h in 1:(nal+nbe+nga)) dY1[jb+1,,j,h] = colSums(dVv[ind,,h])
			}	
		}
		for(h in 1:(npar)){
			count = 0
			for(c in 1:k) for(j in 1:r){
				count = count+1
				ind = (count-1)*mb+(1:mb)
				Cor[h,ind] = t(Am[[j]])%*%(dY1[1:(b[j]+1),c,j,h]-sum(dY1[1:(b[j]+1),c,j,h])*Psi[1:(b[j]+1),c,j])
			}
		}
		for(h in 1:(npar)){
			out = est_multilogit(dV[,,1,h],XXdis,Xlab,be,Pivdis,fort=fort,ex=TRUE)
			Cor[h,nal+(1:nbe)] = out$sc
		}
		dU = array(U,c(k,k,ns,TT,npar))*dlU
        if(param=="multilogit"){
   	    	rGa = dim(Ga)[1]
		    for(h in 1:k){
	    		for(h1 in 1:npar){
				    UU = NULL
			    	for(t in 2:TT) UU = rbind(UU,t(dU[h,,,t,h1]))
			    	out = est_multilogit(UU,ZZdis[,,,h],Zlab,Ga[,h],PIdis[,,h],fort=fort,ex=TRUE)
			    	ind = nal+nbe+(h-1)*rGa+(1:rGa)
			    	Cor[h1,ind] = out$sc
				}
		    }
		}else if(param=="difflogit"){
   	    	rGa = length(Ga)
			for(h1 in 1:npar){
    		   	Tmp = aperm(dU[,,,2:TT,h1],c(1,3,4,2))
			   	Tmp = matrix(Tmp,ns*k*(TT-1),k)
		       	out = est_multilogit(Tmp,ZZdis,Zlab,Ga,PIdis,fort=fort,ex=TRUE)
			    ind = nal+nbe+(1:rGa)
			  	Cor[h1,ind] = out$sc
			}
   	    }
 # check score and information
		if(check_der){
			lk0 = lk
			out = lk_obs_latent(th,S,R,b,yv,Am,XXdis,Xlab,ZZdis,Zlab,param,fort)
			print(lk-out$lk)
			sc1 = out$sc
       		lth = length(th)
        	scn = rep(0,lth); Fn = matrix(0,lth,lth)
   	    	for(h in 1:lth){
   		    	thh = th; thh[h] = thh[h]+10^-6
           		outh = lk_obs_latent(thh,S,R,b,yv,Am,XXdis,Xlab,ZZdis,Zlab,param,fort=fort)
           		scn[h] = (outh$lk-lk)*10^6
           		Fn[,h] = -(outh$sc-sc)*10^6
   	 		}
	     	print(round(cbind(sc,sc1,scn,sc-scn),4))
   		    print(round(cbind(diag(Fi-Cor),diag(Fn),diag(Fi-Cor-Fn)),4))
   		    browser()
		}
# Information matrix and standard errors
    	Fi = Fi-Cor
		iFi = ginv(Fi)		
		se = sqrt(diag(iFi))
# Divide parameters
      	sepsi = se[1:nal]
        sebe = se[nal+(1:nbe)]
        sega = se[nal+nbe+(1:nga)]
	}	
# Compute number of parameters
    if(r==1) np = k*mb*r else np = k*sum(b)
    np = np+(k-1)*(nc1+1)
    if(param=="multilogit") np = np+(k-1)*(nc2+1)*k else if(param=="difflogit")  np = np+(k-1)*(nc2+k)
  	aic = -2*lk+np*2
  	bic = -2*lk+np*log(n)
#	out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=lkv,J=J,V=V1,th=th,sc=sc)
# local decoding
    Ul = matrix(0,ns,TT)
    for(i in 1:ns) for(t in 1:TT){
    	Ul[i,t] = which.max(V[i,,t])
    }
    if(all(yv==1)) V1=V

    if(out_se){
    	if(r==1){
	    	psi = as.vector(aperm(Psi,c(1,3,2)))
  	  		dPsi = diag(psi)%*%matrix(aperm(dlPsi,c(1,3,2,4)),c((b+1)*k,nal)) 
    		sePsi = sqrt(diag(dPsi%*%iFi[1:nal,1:nal]%*%t(dPsi)))
    		sePsi = aperm(array(sePsi,c(b+1,r,c)),c(1,3,2))
    		dimnames(sePsi)=list(category=0:b,state=1:k)
    	}else{
    		sePsi = array(NA,dim(Psi))
    		for(j in 1:r) sePsi[1:(b[j]+1)]=0
    		ind = 0
    		for(c in 1:k) for(j in 1:r){
				Tmp = iFi[ind+(1:b[j]),ind+(1:b[j])]
    			ind = ind+b[j]
    			psi = Psi[1:(b[j]+1),c,j]
    			Tmp1 = matrix((diag(psi)-psi%o%psi)[,-1],b[j]+1,b[j])
    			sePsi[1:(b[j]+1),c,j] = sqrt(diag(Tmp1%*%Tmp%*%t(Tmp1)))
    			dimnames(sePsi)=list(category=0:mb,state=1:k,item=1:r)
    		}
    	}
    }
	Be = matrix(be,nc1+1,k-1)
	if (is.null(nameBe)){
		nameBe = c("intercept",paste("X1",1:nc1,sep=""))
	}else{
		nameBe = c("intercept",nameBe)
	}	
	dimnames(Be) = list(nameBe,logit=2:k)
	if(out_se) {seBe = matrix(sebe,nc1+1,k-1); dimnames(seBe) = list(nameBe,logit=2:k)}
	if(param=="multilogit"){
		if(is.null(nameGa)){
			nameGa = c("intercept", paste("X2",1:nc2,sep=""))
		}else{
			nameGa = c("intercept",nameGa)
		}
		if(k>2) {
			Ga = array(as.vector(Ga),c(nc2+1,k-1,k))
			dimnames(Ga) = list(nameGa,logit=2:k,logit=1:k)
		}else if(k==2){ 
			dimnames(Ga) = 	list(nameGa,logit=1:k)
		}
		if(out_se){
			if(k==2){
				seGa = matrix(sega,nc2+1,2)
				dimnames(seGa) = list(nameGa,logit=1:k)
			}else if(k>2){
				seGa = array(as.vector(sega),c(nc2+1,k-1,k))
				dimnames(seGa) = list(nameGa,logit=2:k,logit=1:k)
			}
		}
	}else if(param=="difflogit"){
		Ga0 = Ga
		Ga = vector("list",2)
		seGa = vector("list",2)
		Ga[[1]] = t(matrix(Ga0[1:(k*(k-1))],k-1,k))
		Ga[[2]] = matrix(Ga0[(k*(k-1))+(1:((k-1)*nc2))],nc2,k-1)
		if(is.null(nameGa)){
			nameGa2 = paste("X2",1:nc2,sep="")
		}else{
			nameGa2 = nameGa
		}
		if (k==2) {
			dimnames(Ga[[1]]) = list(intercept=1:k,logit=k)
			dimnames(Ga[[2]])=list(nameGa2,logit=k)
		} else if (k>2){
			dimnames(Ga[[1]]) = list(intercept=1:k,logit=2:k)
			dimnames(Ga[[2]])=list(nameGa2,logit=2:k)
		}
		if(out_se){
			seGa[[1]] = t(matrix(sega[1:(k*(k-1))],k-1,k))
			seGa[[2]] = matrix(sega[(k*(k-1))+(1:((k-1)*nc2))],nc2,k-1)
			if(k==2){
				dimnames(seGa[[1]]) = list(intercept=1:k,logit=k)
				dimnames(seGa[[2]])=list(nameGa2,logit=k)				
			}else if (k>2){
				dimnames(seGa[[1]]) = list(intercept=1:k,logit=2:k)
				dimnames(seGa[[2]])=list(nameGa2,logit=2:k)
			}
		}
	}
	# adjust output
	lk = as.vector(lk)
	if(output){
		dimnames(Piv)=list(subject=1:n,state=1:k)
		dimnames(PI)=list(state=1:k,state=1:k,subject=1:n,time=1:TT)
	}
	if(r==1) dimnames(Psi) = list(category=0:b,state=1:k,item=1) else 		dimnames(Psi)=list(category=0:mb,state=1:k,item=1:r)
	out = list(lk=lk,Be=Be,Ga=Ga,Psi=Psi,np=np,aic=aic,bic=bic,lkv=lkv,
	           call=match.call(),param=param)
	if(out_se){
		out$sePsi = sePsi
		out$seBe = seBe
		out$seGa = seGa
	} 
	# final output
    if(output){
    	out$V1 = V1
    	out$PI = PI
    	out$Piv = Piv
    	out$Ul = Ul
    }
    #cat(" |-------------|-------------|-------------|-------------|\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    class(out)="LMlatent"
	return(out)
}
