cquad_basic <-
function(id, yv, X=NULL, be=NULL, w = rep(1,n), dyn=FALSE){

# SIMPLIFIED VERSION WITH OF THE ECONOMETRIC MODEL
#
# Fit a Conditional Logit model using CML (with time varying number of observation)
# the first column of data is the increasing number of the unit with aggregated data
# Y       : vector of binary response variable of dimension rc x 1 (r=units, c=times)
# X       : matrix of covariates of dimension rc x k (k=number of covariates)
# beta    : vector of first guess of the regression parameters of dimension k
# w       : vector of weights for each subject
# betahat : vector of estimated parameters of dimension k
# se      : vector of estimated standard errors of dimension k
# lk      : log-likelihood value at convergence
# dyn     : TRUE to include a dynamic component 

# preliminaries
	if(is.null(X) & !dyn) stop("nothing to estimate")
	pid = id
	r = length(pid)
	label = unique(pid)
	n = length(label)
	if(is.null(X)) k=0 else{X = as.matrix(X); k = ncol(X)}
	c = max(yv)+1 # number of categories
	if(k>0) Xv = X
# chech variabiliy of the covariates
	if(k>0) for(j in 1:k){
		xj = X[,j]
		flag = TRUE
		for(il in label){
			xv = xj[which(pid==il)] 
			if(max(xv)-min(xv)>0) flag = FALSE
		}
		if(flag) stop("at least one covariate without variability within unit")
	}
# variable names	
	varnames = NULL
	if(k>0){
		if(is.null(colnames(X))) for(j in 1:k) varnames = c(varnames,paste("X",j)) 
		else varnames = colnames(X)
	}
	if(dyn) varnames = c(varnames,"y_lag")
#  starting values
	if(dyn)  be = rep(0,k+1) else  be = rep(0,k)
# check for balanced data
	Tv = rep(0,n)
	ind = id
	for(i in 1:n) Tv[i] = sum(ind==label[i])
	TT0 = max(Tv)
	if(dyn) Tv = Tv-1
	TT = max(Tv)
	balanced = all(Tv==TT)
  	if(balanced){
  		Y = t(matrix(yv,TT0,n))
  		if(k>0) XX = array(Xv,c(TT0,n,k)) 
  		ZZ = sq(TT); sZZ = rowSums(ZZ)
  	}
  	if(balanced) cat("Balanced panel data\n") else cat("Unbalanced panel data\n") 
# iterate until convergence    
	Sc = matrix(0,n,1)
	it = 0; lk = -Inf; lk0 = -Inf
	if(dyn) zero1 = c(rep(0,k),1)
	cat(" |--------------|--------------|--------------|\n")
	cat(" |   iteration  |      lk      |    lk-lko    |\n")
	cat(" |--------------|--------------|--------------|\n")
	while(abs(lk-lk0)>10^-6 | it==0){
		it = it+1; lk0 = lk
		if(dyn) scv = matrix(0,n,k+1) else scv = matrix(0,n,k)
		lk = 0; J = 0
		for(cut_point in 1:(c-1)){
			yd = 1*(yv>(cut_point-1))
			if(balanced) Yd = 1*(Y>(cut_point-1))
			for(i in 1:n){
				if(Tv[i]>1){
					il = label[i]
					y_i = yd[pid==il]
					if(dyn){y_i0 = y_i[1]; y_i = y_i[-1]}
					sui = sum(y_i)
					if(sui>0 & sui<Tv[i]){
						if(balanced) Z = ZZ[sZZ==sui,]
						else Z = sq(Tv[i],sui)
						if(k>0) if(balanced) x_i = as.matrix(XX[,i,]) else x_i = as.matrix(Xv[pid==il,])
						if(dyn){
							if(k>0) x_i = rbind(cbind(x_i[-1,],0),zero1)
			   				else x_i = as.matrix(c(rep(0,Tv[i]),1))
    						if(Tv[i]==2) Z = cbind(Z,y_i0*Z[,1]+Z[,1]*Z[,2])
    						else Z = cbind(Z,y_i0*Z[,1]+rowSums(Z[,1:Tv[i]-1]*Z[,2:Tv[i]]))
    						xb = x_i%*%be
    						den = exp(Z%*%xb)
    						sden = sum(den)
    						y_i = c(y_i,y_i0*y_i[1]+sum(y_i[1:Tv[i]-1]*y_i[2:Tv[i]]))
    						pc_i = exp(y_i%*%xb)/sden
    					}else{
			    			xb = x_i%*%be
    						den = exp(Z%*%xb)
    						sden = sum(den)
    						pc_i = exp(y_i%*%xb)/sden
    					}
    					Zt = t(Z)
    					lk = lk+w[i]*log(pc_i)
    					pp_i = as.vector(den/sden)
    					e_i = Zt%*%pp_i
    					scv[i,] = scv[i,]+w[i]*(t(y_i-e_i)%*%x_i)
    					V_i = Zt%*%diag(pp_i)%*%Z-e_i%*%t(e_i)
    					J = J-w[i]*(t(x_i)%*%V_i%*%x_i)
    				}
    			}
    		}
    	}
    	sc = colSums(scv)
    	iJ = try(solve(J),silent=TRUE)
    	if(inherits(iJ,"try-error")){
    		iJ = ginv(J)
    		warning("Inversion problems of the information matrix")
    	}
    	be = be-iJ%*%sc
    	cat("",sprintf("%12g", c(it,lk,lk-lk0)), "\n", sep = " | ")
    }
	cat(" |--------------|--------------|--------------|\n")
   	be = as.vector(be)
   	Va = -iJ
   	Var = iJ%*%(t(scv)%*%scv)%*%iJ
   	se = sqrt(abs(diag(Va)))
   	if(any(diag(Va)<0)) warning("Negative elements in asymptotic variance-covariance matrix")
   	ser = sqrt(abs(diag(Var)))
   	if(any(diag(Var)<0)) warning("Negative elements in robust variance-covariance matrix")
   	lk = as.vector(lk)
   	names(be) = varnames
   	colnames(Va) = rownames(Va) = varnames   	
   	colnames(scv) = varnames
   	rownames(J) = colnames(J) = varnames
   	names(se) = varnames   	
   	out = list(formula=formula,lk=lk,coefficients=be,vcov=Va,scv=scv,J=J,se=se,ser=ser,Tv=Tv,call=match.call())
	class(out) = c("cquad","panelmodel")
   	return(out)

}
