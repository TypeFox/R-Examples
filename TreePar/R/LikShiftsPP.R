LikShiftsPP<-function(x,t,lambda,mu,sampling=1,survival=1,root=1,n=0){
	l<-lambda
	if (n==0){
	out<- -log((Ffuncshift(x[1],t,l,mu,sampling)^(root+1)))
	for (i in 2:length(x)){
		out<-out+log(Fderifuncshift(x[i],t,l,mu,sampling))
	}
	out <- out+(1-survival)*(root+1)*log((1-pnshift(0,x[1],t,l,mu,sampling)))}
	else {out<-HelpShiftsPPn(x,t,l,mu,sampling,survival,root)}
	out<- -out
	out
}