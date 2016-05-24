LikShifts <-
function(x,t,lambda,mu,sampling,posdiv=FALSE,survival=1,groups=0) {
	l<-lambda
	res<- -10^12
	boundary<-0
	for (i in 1:length(l)){
		if (l[i]==mu[i]|| mu[i]<0 || l[i]<0.0001 || l[i]>100 || (abs(l[i]-mu[i])<0.0001) ){boundary<-1}
		if (posdiv==TRUE && (l[i]-mu[i])<0.0001 ) {boundary<-1}
		}
	for (i in 1:length(sampling)){
		if (sampling[i]>1 || sampling[i]<=0 ){boundary<-1}
		}
	if (boundary==0) {
	x<-sort(x)
	n<-lineages(x,t)
	mrca<-x[length(x)]
	res<- n[1]*log(sampling[1])  + n[1]* log((l[1]-mu[1])^2) + 2 * log(g(mrca,t,l,mu,sampling))
	if (survival == 1){
		res<-res -2*log(1-q2(inter(mrca,t),mrca,t,l,mu,sampling  ) )}
	for (j in 1:(length(x)-1)) {
		res <- res +log(2*l[inter(x[j],t)]) + log(g(x[j],t,l,mu,sampling))
		}
	if (inter(mrca,t)>1) {
		for (j in 2:inter(mrca,t)) {
			res <- res + n[j] * log(sampling[j] *(l[j]-mu[j])^2 * g((t[j]),t,l,mu,sampling) )
		}
	}
	
	if (length(groups)>1) {
		for (i in 1:length(groups[,1])){
		age<-groups[i,1]
		size<-groups[i,2]
		res <- res + log(pnshift(size,age,t,l,mu))-log(pnshift(1,age,t,l,mu))
		}	
	}
	res<-res-(length(x)-1)*log(2)
	}
	-res
	}

