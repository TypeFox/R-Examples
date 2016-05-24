emis4ps.MSAR <-
function(y,data,theta,covar=NULL) {
	
	d <- attributes(theta)$NbComp
	M <- attributes(theta)$NbRegimes
	order <- attributes(theta)$order
	label <- attributes(theta)$label
	data = as.matrix(data)
	A=NULL
	if(d==1){
		if(order>0){
			for (j in 1:M) {
				A[[j]] = sapply(theta$A[j,],as.numeric)
				   }
		}else{ A[[1]] <- array(0,c(1,1,M)) }
		Sigma=matrix(theta$sigma,M,1)
		A0=matrix(theta$A0,M,1)
	}
	else{
		if(order>0){A=theta$A
		}else{A <- matrix(0,M,d)}	# Inutile? 
		Sigma=theta$sigma
		A0=theta$A0
	}

	T <- dim(as.matrix(data))[1]
	prob <- array(0,c(length(y),T-order,M))
	for (j in 1:M) {
	  	A0.emis = array(0,c(d,T-order))
		if (order>0) {
			for (o in (1:order)) {
				A0.emis = A0.emis + matrix(A[[j]][[o]],d,d) %*% t(data[((order+1):T)-o,])
			}
		}
		A0.emis = A0.emis + A0[j,];

		if (substr(label,2,2)=='N') { 
			nh.emissions <- attributes(theta)$nh.emissions
			par.emis = theta$par.emis
			ncov = dim(covar)[2]
			if (is.null(covar)) {ncov=1}
			for (i in 1:d) {
				femis = nh.emissions(matrix(covar[(order+1):T,],T-order,ncov), as.matrix(par.emis[[j]]))
				A0.emis[i,]=A0.emis[i,]+femis[i,]	
			}
		}
		for (iy in 1:length(y)) {
			prob[iy,,j] = pdf.norm(matrix(y[iy],d,length(A0.emis)),matrix(A0.emis,d,length(A0.emis)),matrix(Sigma[[j]],d,d))
		}
   	}
    prob
}
