est_multi_glob <- function (Y, X, model, ind = 1:nrow(Y), be = NULL, Dis = NULL, dis = NULL, 
                            disp=FALSE, only_sc = FALSE, Int = NULL, der_single = FALSE){

# preliminaries                           	
	if(is.null(dis) & !is.null(Dis)) dis = rep(0,nrow(Dis))
    l = nrow(X) + 1
    ncov = ncol(X)
    n = nrow(Y)
    w = rowSums(Y)
    if (model == "g" || model == "l") {
        out = matr_glob(l, model)
        Co = out$Co
        Ma = out$Ma
    }
    G = rbind(matrix(0, 1, l - 1), diag(l - 1))
    nd = max(ind)
    if (is.null(be)) {
        dist = colMeans(Y)
        if (model == "m") 
            eta = log(dist[-1]/dist[1])
        if (model == "g" || model == "l") 
            eta = Co %*% log(Ma %*% dist)
        Xm = 0
        for (h in 1:nd) Xm = Xm + X[, , h] * sum(w[ind == h])
        Xm = matrix(Xm/sum(w), l - 1, ncov)
        be = ginv(t(Xm) %*% Xm) %*% t(Xm) %*% eta
        if(model=="g"){
        	YW = Y*w
        	num = den = 0
   	        for (h in 1:nd){
		        Xh = matrix(X[, , h], l - 1, ncov)    	
	            indh = ind == h
	            if(sum(indh)==1){
	            	dist = YW[indh,]+0.5
	            }else{
	            	dist = colSums(YW[indh,])+0.5
	            }
	            dist = rep(1,length(dist))
	            eta = Co %*% log(Ma %*% dist)
            	num = num+(t(Xh) %*% eta)*sum(w[indh])
            	den = den+(t(Xh)%*%Xh)*sum(w[indh])
#		        be = be+ginv(den)%*%num*sum(w[indh])/sum(w)
	        }
	        be = ginv(den)%*%num
        }
    }
    P = matrix(0, nd, l)
    for (h in 1:nd) {
        Xh = matrix(X[, , h], l - 1, ncov)
        if(is.null(Int)) x0 = 0 else x0 = Int[,h]
        if (model == "m") {
            P[h, ] = exp(G %*% (x0+Xh %*% be))
            P[h, ] = P[h, ]/sum(P[h, ])
        }
        if (model == "g" | model == "l") 
            P[h, ] = inv_glob(x0+Xh %*% be, model)$p
    }
    if(any(P<0)) print("inversion error in P")
    P = pmax(P, 10^-100)
    Pdis = P; P = P[ind,]
    lk = sum(Y * log(P)) 
    lko = lk
    it = 0
    flag = TRUE
    if(disp) print(c(0,lk))
    while ((abs(lk - lko) > 10^-6 | it == 0) & it < 10^4) {
        it = it + 1
        s = 0; FI = 0
        if(der_single) Sc = matrix(0,n,ncov)
        for (h in 1:nd) {
            indh = which(ind == h)
            p = Pdis[h, ]
            Xh = matrix(X[, , h], l - 1, ncov)
            if (model == "m") 
                D = t(G%*%Xh)
            if (model == "g" || model == "l") {
                ve = 1/as.vector(Ma %*% p)
                D = ginv(Co %*% diag(ve) %*% Ma %*% diag(p) %*% G)
                D = t(G%*%D%*%Xh)
            }
            if(der_single){
	        	for (i in indh){
	        		Sc[i,] = D %*% (Y[i, ] - w[i] * p)     
	        		s = s + Sc[i,]
	        	}       	
            }else{
            	if(length(indh)==1) tmp = Y[indh,]
            	else tmp = as.vector(colSums(Y[indh,]))
	            s = s + D %*% (tmp -sum(w[indh]) * p)
	        }
            FI = FI + sum(w[indh]) * (D %*% (diag(p) - p%*%t(p)) %*% t(D))
        }
        if(!only_sc){
		    if (rcond(FI) < 10^-15){
    			if(flag){
			    	flag = FALSE
    				warning("matrix close to singularity in est_multi_glob")
    			}
        		dbe = ginv(FI) %*% s
        	}
        	else dbe = solve(FI) %*% s
        	mdbe = max(abs(dbe))
        	if (mdbe > 0.25) dbe = dbe/mdbe * 0.25
	        if(!is.null(Dis)){
    	    	LL = chol(FI)
				dbe = lsei(A=LL,B=LL%*%(be+dbe),G=Dis,H=dis,verbose=FALSE)$X-be
       		}
        	be0 = be; be = be + dbe
        	P0 = P; P = matrix(0, nd, l)
        	for (h in 1:nd) {
   		        if(is.null(Int)) x0 = 0 else x0 = Int[,h]
	            if (model == "m") 
    	            P[h, ] = exp(G %*% (x0+X[, , h]%*%be))
	        	    P[h, ] = P[h, ]/sum(P[h, ])
            	if (model == "g" || model == "l") 
                	P[h, ] = inv_glob(x0+X[, , h] %*% be, model)$p
        	}
        	P = pmax(P, 10^-100)
	    	Pdis = P; P = P[ind,]
        	lko = lk
        	lk = sum(Y * log(P))
        	if(is.na(lk)){
        		lk = lko; be = be0; P = P0
        		warning("convergence problemsin est_multi_glob")
        	}
        	if(disp) print(c(it,lk,lk-lko))
        }
    }
    if(der_single) out = list(be = be, lk = lk, Pdis = Pdis, P = P, Sc = Sc, sc = s)
    else out = list(be = be, lk = lk, Pdis = Pdis, P = P, sc = s)
}