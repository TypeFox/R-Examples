MillerJack <-
function(x,y=NA){
check<-0
	if((is.null(ncol(x))||ncol(x)==1)&&!is.na(y)){
		check=1
	}
	if(max(dim(x)[2],1)==2){
		y<-x[!is.na(x[,2]),2]
		x<-x[!is.na(x[,1]),1]
		check=1
	}

	if(!check){
		return('Error: invalid form for entered data')
	}
	outp<-list()
	outp$stat.name<-"Miller Q"
	
	outp$m<-length(x)
	outp$n<-length(y)

	A<-rep(0,outp$m)
	B<-rep(0,outp$n)

	Xbar0<-mean(x)
	S0<-log(sum((x-Xbar0)^2/(outp$m-1)))

	Ybar0<-mean(y)
	T0<-log(sum((y-Ybar0)^2/(outp$n-1)))

	for(i in 1:outp$m){
		Xdel<-x[-i]
		A[i]<-(outp$m*S0-(outp$m-1)*log(var(Xdel)))
	}

	Abar<-mean(A)

	V1<-sum((A-Abar)^2/(outp$m*(outp$m-1)))

	for(i in 1:outp$n){
		Ydel<-y[-i]
		B[i]<-(outp$n*T0-(outp$n-1)*log(var(Ydel)))
	}
	Bbar<-mean(B)
	V2<-sum((B-Bbar)^2/(outp$n*(outp$n-1)))

	obs.stat<-(Abar-Bbar)/(V1+V2)^(1/2)
	return(obs.stat)
}
