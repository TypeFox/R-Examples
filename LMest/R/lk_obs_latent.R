lk_obs_latent <-
function(th,S,R,b,yv,Am,XXdis,Xlab,ZZdis,Zlab,param,fort=TRUE){
	
# preliminaries
   	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	l = max(S)+1
  	if(length(sS)==2) r = 1 else r = sS[3]
  	k = dim(XXdis)[1]
  	nc1 = dim(XXdis)[2]/(k-1)-1
  	if(param=="multilogit"){ 
  		nc2 = dim(ZZdis)[2]/(k-1)-1
  	}else if(param=="difflogit"){
  		nc2 = (dim(ZZdis)[2]-(k*(k-1)))/(k-1)
  	}	
  	mb = max(b)
  	Xndis = max(Xlab)
  	Zndis = max(Zlab)
  	sb = sum(b)
# separate parameters
# conditional response probabilities
  	Psi = array(NA,c(mb+1,k,r))
  	for(j in 1:r) Psi[1:(b[j]+1),,j] = 0 
  	count = 0
    for(c in 1:k) for(j in 1:r){
    	ind = count+(1:b[j])
    	temp = exp(Am[[j]]%*%th[ind])
    	Psi[1:(b[j]+1),c,j] = temp/sum(temp)
    	count = count+b[j]
	}
# parameters on initial probabilities
    ind = count + (1:((1+nc1)*(k-1)))
 	be = th[ind]
   	out = prob_multilogit(XXdis,be,Xlab,fort)
   	Piv = out$P; Pivdis = out$Pdis
   	count = count+((1+nc1)*(k-1))
# parameters on transition probabilities
    if(param=="multilogit"){
    	ind = count+(1:((nc2+1)*(k-1)*k))
        Ga = matrix(th[ind],(nc2+1)*(k-1),k)
		PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,ns,TT))
		for(h in 1:k){
		    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,ns,TT-1))
		}
		count = count+(nc2+1)*(k-1)*k
	}else if(param=="difflogit"){
		ind = count+(1:(k*(k-1)+(k-1)*nc2))
		Ga = matrix(th[ind],k*(k-1)+(k-1)*nc2)
        PI = array(0,c(k,k,ns,TT))
        out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,ns,TT-1))
   		PI = aperm(PI,c(2,1,3,4))
   		count = count+(k*(k-1)+(k-1)*nc2)
	}
# compute corresponding derivatives
    dlPsi = array(NA,c(mb+1,k,r,k*sb))
    for(j in 1:r) dlPsi[1:(b[j]+1),,j,] = 0
    count = 0
	for(c in 1:k) for(j in 1:r){
		ind = count+(1:b[j])
	    temp = pmax(Psi[1:(b[j]+1),c,j],10^-50)
	 	dlPsi[1:(b[j]+1),c,j,ind] = (diag(b[j]+1)-rep(1,b[j]+1)%o%temp)%*%Am[[j]]
	   	count = count+b[j]
	#   	th = c(th,log(temp[-1]/temp[1]))
	}
	dlPiv = array(0,c(ns,k,(1+nc1)*(k-1)))
	for(j in 1:Xndis){
	    temp = pmax(Pivdis[j,],10^-50)
	    Temp = (diag(k)-rep(1,k)%o%temp)%*%XXdis[,,j]
	   	for(i in which(Xlab==j)) dlPiv[i,,] = Temp
    }
  #  th = c(th,be)
    count = 0
    if(param=="multilogit"){
    	dlPI = array(0,c(k,k,ns*TT,(1+nc2)*(k-1)*k))
	    for(h in 1:k){
		 	ind = count+(1:((1+nc2)*(k-1)))
		  	for(j in 1:Zndis){
		   		temp = pmax(PIdis[j,,h],10^-50)
                Temp = (diag(k)-rep(1,k)%o%temp)%*%ZZdis[,,j,h]
                for(i in which(Zlab==j)) dlPI[h,,ns+i,ind] = Temp
		   	}
		   	count = count+((1+nc2)*(k-1))
        #    th = c(th,Ga[,h])		    
		}
        dlPI = array(dlPI,c(k,k,ns,TT,(1+nc2)*(k-1)*k))
	}else if(param=="difflogit"){
		   dlPI = array(0,c(k,k*ns*TT,(k+nc2)*(k-1)))
		   temp0 = rep(1,k); Temp0 = diag(k)		   
		   for(j in 1:Zndis){		   	 
		    	 temp = pmax(PIdis[j,],10^-50)
           	     Temp = (Temp0-temp0%o%temp)%*%ZZdis[,,j]
           	     for(i in which(Zlab==j)) dlPI[,k*ns+i,] = Temp          	 
		    	}
#		    Ga = c(as.vector(t(Ga[[1]])),as.vector(Ga[[2]]))
		  #  th = c(th,Ga)		    
            dlPI = array(dlPI,c(k,k,ns,TT,(k+nc2)*(k-1)))
            dlPI = aperm(dlPI,c(2,1,3,4,5))
	}
# compute log-likelihood
   	out = lk_comp_latent(S,R,yv,Piv,PI,Psi,k,der=TRUE,fort=fort,dlPsi=dlPsi,dlPiv=dlPiv,dlPI=dlPI)
	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv; sc = out$dlk
# output
    out = list(lk=lk,sc=sc,Psi=Psi,be=be,Ga=Ga,Piv=Piv,PI=PI,dlPsi=dlPsi,dlPiv=dlPiv,dlPI=dlPI)
}
