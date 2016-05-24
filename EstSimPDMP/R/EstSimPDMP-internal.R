.Transition <-
function(x,t){
	c<-0;
	while (c==0){
		u<-rnorm(1,mean=x,sd=2+t^2); if (u>0 & u<10){c<-1}
	}
	u
}
.tstar <-
function(x){
	10-x
}
.JumpRate <-
function(x,t){
	1/(1+x)
}
.Transition.DC <-
function(x,t){
	if (x==1){
		if (t<2) {u<-1} else {v<-runif(1); if (v<1/2){u<-2} else {u<-3}}
	}
	else if (x==2){
		v<-runif(1)
		if (v<2/3){u<-1} else {u<-2}
	}
	else if (x==3){
		v<-runif(1); w<-runif(1)
		if ((v<1/2) & (w<1/2)) {u<-3} else {if (v<1/3){u<-1} else {u<-2}}
	}
	u
}
.tstar.DC <-
function(x){
	if (x==1){z<-5} else if (x==2) {z<-6} else if (x==3) {z<-3.5}
	z
}
.JumpRate.DC <-
function(x,t){
	if (x==1){z<-t} else if (x==2) {z<-1} else if (x==3) {z<-t^2}
	z
}
.InvYn <-
function(dat,t){
b=dat[dat>=t]
if (length(b)>0){
	z<-1/length(b)
} else {z<-0}
z
}
.ker <-
function(x){
	z<-c()
	for (k in x){
		if (abs(k)<1){
			z<-c(z,(3/4)*(1-k^2))
		}
		else {z<-c(z,0)}
	}
	z
}
.Tri1 <-
function(dat , x){
	N<-length(dat[1,])-1; m<-length(dat[,1]);
	
	A<- ( dat[,1:N]==matrix(x,nrow=m,ncol=N,byrow=TRUE))
	
	B<-A %*% as.matrix( rep(1,N)); C<-which(B==N);
	
	dat[C,N+1]
}
.Tri2 <-
function(dat , x , y){
	
	N<-length(dat[1,])-1; m<-length(dat[,1]);
	
	A<- ( dat[1:(m-1),1:N]==matrix(x,nrow=(m-1),ncol=N,byrow=TRUE))
	B<- ( dat[2:m,1:N]==matrix(y,nrow=(m-1),ncol=N,byrow=TRUE))	
	
	A2<-A %*% as.matrix( rep(1,N));
	B2<-B %*% as.matrix( rep(1,N));
	C<-which( (A2==N) & (B2==N) );
	
	dat[C,N+1]
}
.CondSurv <-
function(dat , x , y , t){
	A<-.Tri1(dat , x); B<-.Tri2(dat , x , y);
	if (length(A)>0){
		z<-length(B[B>t])/length(A)
	} else {z<-0}
	z
}
.InvYnBis <-
function(dat,t){
	b<-dat[dat<=t]
	z<-c()
	if (length(b)>0){
		for (k in 1:length(b)){
			z<-c(z,.InvYn(dat,b[k]))
		}
	} else {z<-c()}
	z
}
