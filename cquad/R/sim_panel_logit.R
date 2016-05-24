sim_panel_logit <-
function(id,al,X=NULL,eta,dyn=FALSE){
	
# simulate data from the static/dynamic logit model given a vector of ordered indices id,
# a vector of individual effects al of dimension n, a matrix of covariates X,
# and a vector of parameters eta including the parameter for state depedence
# dyn = TRUE to simulate from the dynamic logit model, case in which the first observation
# is simulated from the model with state dependence parameter 0

# preliminaries
	n = max(id); nt = length(id)
	if(is.null(X)) nc = 0 
	else if(is.vector(X)){nc = 1; X = as.vector(X)} else nc = ncol(X)
	if(dyn){
		if(nc!=length(eta)-1) stop("size of eta non-conform with that of X")
	}else{
		if(nc!=(length(eta))) stop("size of eta non-conform with that of X")		
	}
# simulate data
	yv = pv = rep(0,nt)
	if(!dyn){
		for(i in 1:n){
			ind = which(id==i); Ti = length(ind)
			if(nc==0) p = exp(al[i])
			else if(nc==1) p = exp(al[i]+X[ind]*eta)
			else p = exp(al[i]+X[ind,]%*%eta)
			p = p/(1+p)
			yv[ind] = 1*(runif(Ti)<p)
			pv[ind] = p	
		}
	}else{
		for(i in 1:n){
			ind = which(id==i); Ti = length(ind)
			y = p = rep(0,Ti)
			if(nc==0) p[1] = exp(al[i])
			else if(nc==1) p[1] = exp(al[i]+X[ind[1]]%*%eta[1]) 
			else p[1] = exp(al[i]+X[ind[1],]%*%eta[1:nc])
			p[1] = p[1]/(1+p[1])
			y[1] = 1*(runif(1)<p[1])
			if(Ti>1) for(t in 2:Ti){
				if(nc==0) p[t] = exp(al[i]+y[t-1]*eta[nc+1])
				else if(nc==1) p[t] = exp(al[i]+X[ind[t]]%*%eta[1]+y[t-1]*eta[2]) 
				else p[t] = exp(al[i]+X[ind[t],]%*%eta[1:nc]+y[t-1]*eta[nc+1]) 
				p[t] = p[t]/(1+p[t])
				y[t] = 1*(runif(1)<p[t])		
			}
			yv[ind] = y; pv[ind] = p	
		}		
	}		
# final output
	out = list(yv=yv,pv=pv)
	return(out)

}
