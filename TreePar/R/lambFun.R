lambFun<-function(k,N,lamb0,model){
	if (model== 0){f <- lamb0}
	else if (model== -1){
		if (k>N) {f<-0}
		else {
		f <- lamb0*(1-k/N)}}
	else {	f <- lamb0*k^(-model)}
	f
}
