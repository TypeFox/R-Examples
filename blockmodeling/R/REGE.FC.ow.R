"REGE.FC.ow" <-
function(M,E=1,iter=3,until.change=TRUE,use.diag=TRUE,normE=FALSE){
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
				den<-0
				for(ir in 1:nr){
					for(k in 1:n){
						if(M[i,k,ir]>0) {
							mins<-Eall[k,,it]*pmin(M[i,k,ir],M[j,,ir])
							num<-num+max(mins)
							den<-den+min(pmax(M[i,k,ir],M[j,which(mins==max(mins)),ir]))
							#cat("M[i,k]: ","i = ",i, ", j = ",j,", k = ",k,", num = ", num, ", den = ", den,"\n")
						}

						if(M[k,i,ir]>0) {
							mins<-Eall[k,,it]*pmin(M[k,i,ir],M[,j,ir])
							num<-num+max(mins)
							den<-den+min(pmax(M[k,i,ir],M[which(mins==max(mins)),j,ir]))
							#cat("M[k,i]: ","i = ",i, ", j = ",j,", k = ",k,", num = ", num, ", den = ", den,"\n")
						}

						if(M[j,k,ir]>0) {
							mins<-Eall[k,,it]*pmin(M[j,k,ir],M[i,,ir])
							num<-num+max(mins)
							den<-den+min(pmax(M[j,k,ir],M[i,which(mins==max(mins)),ir]))
							#cat("M[j,k]: ","i = ",i, ", j = ",j,", k = ",k,", num = ", num, ", den = ", den,"\n")
						}

						if(M[k,j,ir]>0) {
							mins<-Eall[k,,it]*pmin(M[k,j,ir],M[,i,ir])
							num<-num+max(mins)
							den<-den+min(pmax(M[k,j,ir],M[which(mins==max(mins)),i,ir]))
							#cat("M[k,j]: ","i = ",i, ", j = ",j,", k = ",k,", num = ", num, ", den = ", den,"\n")
						}
					}
				}
				if(den!=0) {
					Eall[j,i,it+1]<-Eall[i,j,it+1]<-num/den
				} else Eall[j,i,it+1]<-Eall[i,j,it+1]<-1
					diag(Eall[,,it+1])<-1
				}
			}
			if(normE){
				diag(Eall[,,it+1])<-0
				Eall[,,it+1]<-Eall[,,it+1]/sqrt(outer(apply(Eall[,,it+1],1,sum), apply(Eall[,,it+1],2,sum)))
				diag(Eall[,,it+1])<-max(Eall[,,it+1])
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
