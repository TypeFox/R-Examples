"condno" <-
function(W, type){

W <- as.matrix(W)
Winv <- Rmatsolve(W)
if (type=="l1"){
	n1<-sum(abs(W))
	n2<-sum(abs(Winv))
}

if(type=="1" | type=="o" | type=="O"){
	n1<-max(abs(apply(W,2,sum)))
	n2<-max(abs(apply(Winv,2,sum)))
}
if(type=="i" | type=="I"){
	n1<-max(abs(apply(W,1,sum)))
	n2<-max(abs(apply(Winv,1,sum)))
}
if(type=="f" | type=="F"){
	q<-as.vector(W)
	r<-as.vector(Winv)
	n1<-sum(abs(q)^2)
	n2<-sum(abs(r)^2)
}
if(type=="m" | type=="M"){
	n1<-max(abs(W))
	n2<-max(abs(Winv))
}

condno <- n1 * n2
condno
}
