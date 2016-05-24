mnntsmarginal <- function(cestimatesarray,M,component,theta){

cparameters <- matrix(0,nrow=M[component]+1,ncol=M[component]+1)

for (k in 0:M[component]){
	for (m in 0:M[component]){
		cparameters[k+1,m+1] <- ((2*pi)^(length(M)-1))*sum(cestimatesarray[cestimatesarray[,component]==k,length(M)+1]*Conj(cestimatesarray[cestimatesarray[,component]==m,length(M)+1]))
	}
}

res<-0
for (k in 0:M[component]){
	for (m in 0:M[component]){
		res <- res + cparameters[k+1,m+1]*exp((0+1i)*(k - m)*theta)
	}
}
res<-Re(res)
return(res)

}
