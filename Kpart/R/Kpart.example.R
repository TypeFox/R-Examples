Kpart.example <-
function(data=NULL,K=4)
{

	if(is.null(data)) {
		p = paste(system.file("extdata", package="Kpart"), "/example_data.txt", sep="")
		data = read.table(p, header=T)
	}

	#variable assignment

	x=data[,2]
	x.2=x^2
	x.3=x^3
	y=data[,1]
	L=floor(length(y)/K)

	#Min/Max algorithm to find absolute maximum deviate in the kth partition

	m=matrix(nrow=K,ncol=1)
	for (k in 1:K) m[k]=abs(max(y[(L*(k-1)+1):(k*L)])-mean(y[(L*(k-1)+1):(k*L)]))

	W=matrix(nrow=L, ncol=K)
	for (l in 1:L) for (k in 1:K) W[l,k]=ifelse(m[k] != abs(y[(L*(k-1)+l)]-mean(y[(L*(k-1)+1):(k*L)])), 0, x[(L*(k-1)+l)])
	pops=colSums(W)

	#matrix of potential spline knots

	TT=matrix(nrow=length(x), ncol=K)
	for (i in 1:length(x)) for (k in 1:K) TT[i,k]=(ifelse(0 > (x[i]-pops[k])^3 , 0 ,(x[i]-pops[k])^3))
	bigT=data.frame(TT)
	for(i in 1:length(pops))
	colnames(bigT)[i]=paste("X",as.character(pops[i]),sep="")

	#model selection process

	X=cbind(x,x.2,x.3,bigT)
	r=regsubsets(X,y,,method="exhaustive",nbest=1,nmax=(K+4))
	a=(summary(r)[c("bic")])
	AA=data.frame(a)
	V=cbind(1:length(t(AA)),AA)
	aa=matrix(nrow=length(t(AA)),ncol=1)
	for (i in 1:length(t(AA))) aa[i,1]=ifelse(min(V[,2])==V[i,2],V[i,1],0)
	aaa=data.frame(aa)
	aaaa=colSums(aaa)
	c=coef(r,id=aaaa)
	cc=data.frame(c)
	C=matrix(nrow=2,ncol=length(c))
	C[2,]=cc[,1]
	C[1,]=names(c)

	#matrix of spline knots

	t=matrix(nrow=length(C[1,]),ncol=length(names(bigT)))
	for (i in 1:length(C[1,])) 
	for (j in 1:length(names(bigT))) 
	t[i,j]=(ifelse(C[1,i]==(names(bigT)[j]),(names(bigT)[j]),NA))
	tt=na.omit(c(t))

	P=matrix(nrow=K, ncol=1)
	for (i in 1:length(tt)) P[i]=ifelse(tt[i]=="NA", NA, paste("bigT",tt[i], sep="$"))
	PP=na.omit(as.character(P))

	#matrix of x,x^2, and x^3 terms 

	qq=matrix(nrow=3,ncol=3)
	for (i in 2:4) qq[(i-1),1]=ifelse(C[1,i]=="x","x",NA)
	for (i in 2:4) qq[(i-1),2]=ifelse(C[1,i]=="x.2","x.2",NA)
	for (i in 2:4) qq[(i-1),3]=ifelse(C[1,i]=="x.3","x.3",NA)
	QQ=na.omit(as.character(qq))
	Q=c(QQ)

	EE=c(Q,PP)
	E=paste(EE,collapse="+")

	#paste of terms to be estimated by lm()

	FF=as.formula(paste("y~ ", E))
	f=lm(FF)

	#plot of x and y axis and fitted values with vertical lines at knots 

	plot(x,y,xlab="your_x",ylab="your_y")
	lines(x,f$fitted.values, col="red")

	lineyears=matrix(nrow=length(PP),ncol=length(colSums(W)))
	for(i in 1:length(PP))
	for(j in 1:length(colSums(W)))
	lineyears[i,j]=ifelse(data.frame(strsplit(PP[i],"X"))[2,1]==colSums(W)[j],colSums(W)[j],0)

	spline.knots=colSums(lineyears)[which(colSums(lineyears)!=0)]

	par(new=TRUE)
	for (i in 1:length(spline.knots)) 
	abline(v=spline.knots[i], lty=2)
	#prints spline knots used in model

	print(spline.knots)
	#returns linear model which summary(lm(.)) can be used

	return(f)

}
