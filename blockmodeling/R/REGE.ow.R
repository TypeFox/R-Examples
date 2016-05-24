"REGE.ow" <-
function(M,E=1,iter=3,until.change=TRUE,use.diag=TRUE){
	n<-dim(M)[1]
	if(n!=dim(M)[2]) stop("M must be a 1-mode matrix")
	if(length(dim(M))==2)M<-array(M,dim=c(n,n,1))
	nr<-dim(M)[3]
	if(!use.diag){for(ir in 1:nr) diag(M[,,ir])<-0}

	Eall<-array(NA,dim=c(n,n,iter+1)) #An array of 'iter' similiaritie matrices
	Eall[,,1]<-E
	for(it in 1:iter){
		for(i in 2:n){
			for(j in 1:(i-1)){
				num<-0
				for(ir in 1:nr){
					for(k in 1:n){
						if(M[i,k,ir]>0) {
							num<-num+max(Eall[k,,it]*pmin(M[i,k,ir],M[j,,ir]))
						}

						if(M[k,i,ir]>0) {
							num<-num+max(Eall[k,,it]*pmin(M[k,i,ir],M[,j,ir]))
						}

						if(M[j,k,ir]>0) {
							num<-num+max(Eall[k,,it]*pmin(M[j,k,ir],M[i,,ir]))
						}					

						if(M[k,j,ir]>0) {
							num<-num+max(Eall[k,,it]*pmin(M[k,j,ir],M[,i,ir]))
						}

					}
				}
				den<-sum(M[c(i,j),,])+sum(M[,c(i,j),])
				if(den!=0) {
					Eall[j,i,it+1]<-Eall[i,j,it+1]<-num/den
				} else Eall[j,i,it+1]<-Eall[i,j,it+1]<-1
				diag(Eall[,,it+1])<-1
			}
		}
		
		if(until.change & all(Eall[,,it]==Eall[,,it+1])){
			Eall<-Eall[,,1:(it+1)]
			break
		}
	}
	itnames<-0:(it)
	itnames[1]<-"initial"
	itnames[it+1]<-"final"
	dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],itnames)
	return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter,use.diag=use.diag))
}

