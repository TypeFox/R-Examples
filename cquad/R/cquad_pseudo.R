cquad_pseudo <-
function(id, yv, X=NULL, be=NULL){

# Fit a Conditional Logit model using CML (with time varying number of observation)
# the first column of data is the increasing number of the unit with aggregated data
# Y       : vector of binary response variable of dimension rc x 1 (r=units, c=times)
# X       : matrix of covariates of dimension rc x k (k=number of covariates)
# beta    : vector of first guess of the regression parameters of dimension k
# w       : vector of weights for each subject
# betahat : vector of estimated parameters of dimension k
# se      : vector of estimated standard errors of dimension k
# lk      : log-likelihood value at convergence

# preliminaries
	pid = id
	r = length(pid)
	label = unique(pid)
	n = length(label)
	if(is.null(X)) k=0 else{X = as.matrix(X); k = ncol(X)}
	if(k>0) Xv = X
# chech variabiliy of the covariates
	if(k>0) for(j in 1:k){
		flag = TRUE
		for(i in 1:n){
			il = label[i] 
			if(max(X[pid==il,j])-min(X[pid==il,j])>0) flag = FALSE
		}
		if(flag) stop("at least one covariate without variability within unit")
	}
# variable names	
	varnames = NULL
	if(k>0){
		if(is.null(colnames(X))) for(j in 1:k) varnames = c(varnames,paste("X",j,sep="")) 
		else varnames = colnames(X)
	}
	varnames = c(varnames,"y_lag")	
# check for balanced data
	Tv = rep(0,n)
	ind = id
	for(i in 1:n) Tv[i] = sum(ind==label[i])
	TT0 = max(Tv)
	Tv = Tv-1
	TT = max(Tv)
	balanced = all(Tv==TT)
  	if(balanced){
  		Y = t(matrix(yv,TT0,n))
  		if(k>0) XX = array(Xv,c(TT0,n,k)) 
  		ZZ = sq(TT); sZZ = rowSums(ZZ)
  	}
#  starting values
	if(k==0){
		be = 0
		Q = rep(0.5,length(yv))
	}else{
		cat("First step estimation\n")
		out = cquad_basic(id,yv,X)
		be0 = be = out$coefficients; scv0 = out$scv; J0 = out$J
		Q = rep(0,length(yv))
		for(i in 1:n){
			il = label[i]
			y_i = yv[pid==il]
			if(all(y_i==0)) q_i = rep(0,length(y_i))
			else if(all(y_i==1)) q_i = rep(1,length(y_i))
			else{
				if(k==0) x_i = NULL else x_i = as.matrix(Xv[pid==il,])
				int = x_i%*%be
				al = 0
	    		q_i = exp(int); q_i = q_i/(1+q_i)
		    	lk1 = as.vector(y_i%*%log(q_i)+(1-y_i)%*%log(1-q_i)); lk1o = -Inf
		   		while(abs(lk1-lk1o)>10^-10){
	    			lk1o = lk1
	    			dal = sum(y_i-q_i)/sum(q_i*(1-q_i))
	    			mdal = abs(dal)
	    			if(mdal>0.5) dal = dal/mdal*0.5
	   				al = al+dal
					q_i = exp(al+int); q_i = q_i/(1+q_i)
					lk1 = as.vector(y_i%*%log(q_i)+(1-y_i)%*%log(1-q_i))
				}
			}
	   		Q[pid==il] = q_i
		}
		be = c(be,0)
	}
# iterate until convergence    
	Sc = matrix(0,n,1)
	it = 0; lk = -Inf; lk0 = -Inf
	zero1 = c(rep(0,k),1)
	cat("Second step estimation\n")
	cat(" |--------------|--------------|--------------|\n")
	cat(" |   iteration  |      lk      |    lk-lko    |\n")
	cat(" |--------------|--------------|--------------|\n")
	while(abs(lk-lk0)>10^-6 | it==0){
		it = it+1; lk0 = lk
		scv = matrix(0,n,k+1)
		lk = 0; J = 0
		for(i in 1:n){
			if(Tv[i]>1){
				il = label[i]
				y_i = yv[pid==il]; y_i0 = y_i[1]; y_i = y_i[-1]
				q_i = Q[pid==il]; q_i = q_i[-1]
				sui = sum(y_i)
				if(sui>0 & sui<Tv[i]){
					if(balanced) Z = ZZ[sZZ==sui,]
					else Z = sq(Tv[i],sui)
					if(k==0) x_i = NULL else x_i = as.matrix(Xv[pid==il,])
					if(k>0) x_i = rbind(cbind(x_i[-1,],0),zero1)
	   				else x_i = as.matrix(c(rep(0,Tv[i]),1))
   					if(Tv[i]==2) Z = cbind(Z,y_i0*(Z[,1]-q_i[1])+Z[,1]*(Z[,2]-q_i[2]))
   					else Z = cbind(Z,y_i0*(Z[,1]-q_i[1])+
   								 rowSums(Z[,1:Tv[i]-1]*(Z-rep(1,nrow(Z))%o%q_i)[,2:Tv[i]]))
#   								 	   					if(it==1 & i==1) print(Z)

           			xb = x_i%*%be
           			den = exp(Z%*%xb)
           			sden = sum(den)
           			y_i = c(y_i,y_i0*(y_i[1]-q_i[1])+sum(y_i[1:Tv[i]-1]*(y_i-q_i)[2:Tv[i]]))
           			pc_i = exp(y_i%*%xb)/sden
       				Zt = t(Z)
       				lk = lk+log(pc_i)
       				pp_i = as.vector(den/sden)
       				e_i = Zt%*%pp_i
       				scv[i,] = scv[i,]+(t(y_i-e_i)%*%x_i)
       				V_i = Zt%*%diag(pp_i)%*%Z-e_i%*%t(e_i)
       				J = J-(t(x_i)%*%V_i%*%x_i)
       			}
       		}
       	}
		sc = colSums(scv)
		iJ = solve(J)
		dbe = -iJ%*%sc
    	mdbe = max(abs(dbe))
   		if(mdbe>0.5) dbe = dbe/mdbe*0.5
		be = be+dbe
   		cat("",sprintf("%12g", c(it,lk,lk-lk0)), "\n", sep = " | ")
	}
	cat(" |--------------|--------------|--------------|\n")
	be = as.vector(be)
	Va = iJ%*%(t(scv)%*%scv)%*%t(iJ)
	if(k==0){
		Va2 = Va
	}else{
		scv1 = cbind(scv0,scv)
# numerical derivative for mixed second derivative
# ------------------------------------------------	
		J10 = matrix(0,k+1,k)
		for(j in 1:k){
			be1 = be0; be1[j] = be1[j]+10^-6
			Q1 = rep(0,length(yv))
			for(i in 1:n){
				il = label[i]
				y_i = yv[pid==il]
				if(all(y_i==0)) q_i = rep(0,length(y_i))
				else if(all(y_i==1)) q_i = rep(1,length(y_i))
				else{
					if(k==0) x_i = NULL else x_i = as.matrix(Xv[pid==il,])
					int = x_i%*%be1
					al = 0
	   				q_i = exp(int); q_i = q_i/(1+q_i)
	    			lk1 = as.vector(y_i%*%log(q_i)+(1-y_i)%*%log(1-q_i)); lk1o = -Inf
	    			while(abs(lk1-lk1o)>10^-6){
	   					lk1o = lk1
		    			dal = sum(y_i-q_i)/sum(q_i*(1-q_i))
		    			mdal = abs(dal)
	    				if(mdal>0.5) dal = dal/mdal*0.5
	   					al = al+dal
						q_i = exp(al+int); q_i = q_i/(1+q_i)
						lk1 = as.vector(y_i%*%log(q_i)+(1-y_i)%*%log(1-q_i))
					}
				}
				Q1[pid==il] = q_i
			}
			sc1 = rep(0,k+1)
			for(i in 1:n){
				if(Tv[i]>1){
					il = label[i]
					y_i = yv[pid==il]; y_i0 = y_i[1]; y_i = y_i[-1]
					q_i = Q1[pid==il]; q_i = q_i[-1]
					sui = sum(y_i)
					if(sui>0 & sui<Tv[i]){
						if(balanced) Z = ZZ[sZZ==sui,]
						else Z = sq(Tv[i],sui)
						if(k==0) x_i = NULL else x_i = as.matrix(Xv[pid==il,])
						if(k>0) x_i = rbind(cbind(x_i[-1,],0),zero1)
		   				else x_i = as.matrix(c(rep(0,Tv[i]),1))
	   					if(Tv[i]==2) Z = cbind(Z,y_i0*(Z[,1]-q_i[1])+Z[,1]*(Z[,2]-q_i[2]))
	   					else Z = cbind(Z,y_i0*(Z[,1]-q_i[1])+
	   								 rowSums(Z[,1:Tv[i]-1]*(Z-rep(1,nrow(Z))%o%q_i)[,2:Tv[i]]))
	           			xb = x_i%*%be
	           			den = exp(Z%*%xb)
	           			sden = sum(den)
	           			y_i = c(y_i,y_i0*(y_i[1]-q_i[1])+sum(y_i[1:Tv[i]-1]*(y_i-q_i)[2:Tv[i]]))
	           			pc_i = exp(y_i%*%xb)/sden
	       				Zt = t(Z)
	       				pp_i = as.vector(den/sden)
	       				e_i = Zt%*%pp_i
	       				sc1 = sc1+(t(y_i-e_i)%*%x_i)
	       			}
	       		}
	       	}
			J10[,j] = (sc1-sc)*10^6		
		}
# ------------------------------------------------	
		J1 = rbind(cbind(J0,matrix(0,k,k+1)),
				  cbind(J10,J))
		iJ1 = solve(J1)
		Va2 = iJ1%*%(t(scv1)%*%scv1)%*%(iJ1)
		Va2 = Va2[-(1:k),-(1:k)]
	}
	se = sqrt(diag(Va))
	se2 = sqrt(diag(Va2))
   	lk = as.vector(lk)
   	names(be) = varnames
   	colnames(Va) = rownames(Va) = varnames   	
   	colnames(scv) = varnames
   	rownames(J) = colnames(J) = varnames
   	names(se) = varnames   	
   	names(se2) = varnames	
	out = list(formula=formula,lk=lk,coefficients=be,vcov=Va,scv=scv,J=J,se=se,ser=se2,Tv=Tv,call=match.call())
	class(out) = c("cquad","panelmodel")
	return(out)

}
