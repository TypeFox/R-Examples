"REGE" <-
function(M,E=1,iter=3,until.change=TRUE,use.diag=TRUE){
	n<-dim(M)[1]
	if(n!=dim(M)[2]) stop("M must be a 1-mode matrix")
	if(!use.diag)diag(M)<-0
	Eall<-array(NA,dim=c(n,n,iter+1)) #An array of 'iter' similiaritie matrices
	Eall[,,1]<-E
	Match<-array(NA,dim=rep(n,4))
	for(i in 2:n){
		for(j in 1:(i-1)){
			for(k in 1:n){
				for(m in 1:n){
					Match[i,j,k,m]<-min(M[i,k],M[j,m]) + min(M[k,i],M[m,j])
					Match[j,i,k,m] <- min(M[j,k],M[i,m]) + min(M[k,j],M[m,i])#/max(1,(max(M[i,k],M[j,m]) + max(M[k,i],M[m,j])+max(M[j,k],M[i,m]) + max(M[k,j],M[m,i])))
				}
			}	
		}
	}

	for(it in 1:iter){
		for(i in 2:n){
			for(j in 1:(i-1)){
				num<-0
				for(k in 1:n){
					#sim<-max(Eall[k,,it]*Match[i,j,k,])
					num<-num+max(Eall[k,,it]*Match[i,j,k,])+max(Eall[k,,it]*Match[j,i,k,])
					#if(i==2&j==1)cat("num = ", num,", den = ",den,", k = ",k,", Maxm1 = ",Maxm1,", ms1 = ",ms1,", Maxm2 = ",Maxm2,", ms2 = ",ms2,"\n")
				}
				#cat("iter=",it,", i=",i,", j=",j,", num=",num,", den=", den,"\n")
				den<-sum(M[c(i,j),])+sum(M[,c(i,j)])
				if(den!=0) {
					Eall[j,i,it+1]<-Eall[i,j,it+1]<-num/den
				} else Eall[j,i,it+1]<-Eall[i,j,it+1]<-1
			}
		}
		diag(Eall[,,it+1])<-1
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

