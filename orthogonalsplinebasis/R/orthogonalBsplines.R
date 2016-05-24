setClass("SplineBasis",representation(knots="numeric", order="integer", Matrices="array"))
setClass("OrthogonalSplineBasis", representation("SplineBasis", transformation="matrix"))

SplineBasis<-function(knots, order=4, keep.duplicates=FALSE) { 
	#degree<-order-1
	n<-length(knots)
	if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
		knots<-unique(knots[order:(n-order+1)])
		knots<-knots[c(1,1,1,seq(length(knots)),length(knots),length(knots),length(knots))]
	}
	q<-n-order
	M<-FindSplineMatrices(knots,order)
	M2<-array(dim=c(order,q,n-2*order+1))
	for(i in order:(n-order)) {  #Identifying interior intervals
		M2[,,i-order+1]<-M[,,i]%*%I_Sel(i,order,n)
	}
	new("SplineBasis", knots=knots, order=as.integer(order), Matrices=M2)
}
GramMatrix<-function(object){
	M<-object@Matrices
	k<-object@order
	Delta<-1/Hankel(1:(2*k-1),k,k)
	di<-diff(object@knots[(k):(length(object@knots)-k+1)])
	d<-dim(M)
	s<-matrix(0,d[2],d[2])
	for(i in 1:d[3]) { 
		s<-s+di[i]*t(M[,,i])%*%Delta%*%M[,,i]
	}
	return(s)
}
OrthogonalizeBasis<-function(object,...){
	s<-GramMatrix(object)
	L<-solve(chol(s))
	M<-object@Matrices
	N<-M
	d<-dim(M)
	for( i in 1:d[3]) { 
		N[,,i]<-M[,,i]%*%L
	}
	return(new("OrthogonalSplineBasis", knots=object@knots, order=as.integer(object@order), Matrices=N, transformation=L))
}
OrthogonalSplineBasis<-function(knots,...)OrthogonalizeBasis(SplineBasis(knots,...))
setGeneric("orthogonalize",function(object,...)standardGeneric("orthogonalize"))
setMethod("orthogonalize",signature("SplineBasis"),OrthogonalizeBasis,valueClass="OrthogonalSplineBasis")
dim.SplineBasis<-function(x)dim(x@Matrices)
setMethod("dim","SplineBasis", dim.SplineBasis)
EvaluateBasis<-function(object,x,...) { 
	stopifnot(is.numeric(x))
	dots<-list(...)
	if(length(x)>1){
		results<-t(sapply(x,EvaluateBasis, object=object))
	} else {
		M<-object@Matrices
		knots<-object@knots
		order<-object@order
		if(x < knots[order] | x>knots[length(knots)-order+1])return(rep(NA,dim(object)[2]))
		if(x==knots[length(knots)-order+1]){ ind <- x<=knots } else { ind <- x<knots }
		if(all(ind)|all(!ind))  {
			if(x==knots[length(knots)-order+1]) 
				return(rep(1,order)%*%matrix(M[,,dim(M)[3]],nrow=order)) 
			else 
				return(rep(0,ncol(object)))
		}
		i<-which(ind)[1]-1
		u<-(x-knots[i])/(knots[i+1]-knots[i])
		U<-u^(0:(order-1))
		return(U%*%matrix(M[,,i-order+1],nrow=order))
	}
} 
setGeneric("evaluate",function(object, x,...)standardGeneric("evaluate"))
setMethod("evaluate",signature("SplineBasis","numeric"),function(object, x, ...)EvaluateBasis(object=object, x=x, ...))
deriv.SplineBasis<-function(object,l=1) {
	dnew<-d<-dim(object)
	DM<-DerivativeMatrix(d[1])
	M<-object@Matrices
	k<-object@order
	knots<- object@knots
	n<-length(knots)
	
	dnew[1]=d[1]-l
	N<-array(0,dim=dnew)
	w<-(diff(knots)[k:(n-k)])^-l
	for( i in 1:d[3]) { 
		N[,,i]<-(w[i]*t(MatrixPower(DM,l))%*%M[,,i])[seq(d[1]-l),]
	}
	
	newknots<-knots[seq((l+1),(n-l))]
	order=object@order-l
	new("SplineBasis",knots=newknots, order=as.integer(order), Matrices=N)
}
setMethod("deriv", signature("SplineBasis"), function(expr,...)deriv.SplineBasis(object=expr,...))
integrate.SplineBasis<-function(object,...){
	dnew<-d<-dim(object)
	M<-object@Matrices
	k<-object@order
	knots<- object@knots
	n<-length(knots)
	
	dnew[1]=d[1]+1
	w<-(diff(knots)[k:(n-k)])
	
	N<-array(0,dim=dnew)
	for(i in 1:d[3]){
		N[,,i]<-rbind( if(i>1)(rep(1,k+1))%*%N[,,i-1] else 0,w[i]*t(diag(1/(1:k)))%*%M[,,i] )
	}
	
	newknots<-knots[c(1,seq(n),n)]
	order=object@order+1L
	new("SplineBasis",knots=newknots, order=as.integer(order), Matrices=N)
}
setGeneric("integrate",function(object,...)standardGeneric("integrate"))
setMethod("integrate",signature("SplineBasis"),integrate.SplineBasis)

print.SplineBasis<-function(object) { 
	cat("Spline Basis\n")
	cat("Order: ",object@order,"\n",
		"Degree: ",object@order-1,"\n",
		"Knots: ",paste(object@knots,collapse=" "),"\n",sep="")
	invisible(object)
}
setMethod("show","SplineBasis",print.SplineBasis)
print.OrthogonalSplineBasis<-function(object) { 
	cat("Orthogonalized ")
	invisible(print.SplineBasis(object))
}
setMethod("show","OrthogonalSplineBasis",print.OrthogonalSplineBasis)
plot.SplineBasis<-function(x,y,xlab=NULL,ylab=NULL,main='Basis Functions', type='l', ...) { 
	dots<-list(...)
	dots$x<-plotdata<-seq(x@knots[x@order],x@knots[length(x@knots)-x@order+1],length=1000)
	dots$y<-evaluate(x,as(plotdata,"numeric"))
	if(is.null(xlab)) xlab<-'' 
	if(is.null(ylab)) ylab<-''
	# if(!hasArg(ylab)) dots$ylab<-"Basis Functions" 
	# if(!hasArg(main)) dots$main<-as.character(substitute(x,as.environment(-1))) 
	# if(!hasArg(type)) dots$type<-'l'
	# browser()
	
	# do.call("matplot",dots)
	matplot(plotdata,evaluate(x,as(plotdata,"numeric")),type='l',xlab=xlab,ylab=ylab,main=main,...)
}
setMethod("plot",signature(x="SplineBasis",y="missing"),plot.SplineBasis)
setMethod("plot",signature(x="SplineBasis",y="vector"),
	function(x,y,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),type='l',...){
	stopifnot(NROW(y)==ncol(x))
	k<-x@knots[(x@order):(length(x@knots)-x@order+1)]
	px<-seq(min(k),max(k),length=100)
	B<-evaluate(x,px)
	plot(px,B%*%y,xlab=xlab,ylab=ylab,type=type,...)
})
setMethod("plot",signature(x="SplineBasis",y="matrix"),
	function(x,y,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),type='l',...){
	stopifnot(NROW(y)==ncol(x))
	k<-x@knots[(x@order):(length(x@knots)-x@order+1)]
	px<-seq(min(k),max(k),length=100)
	B<-evaluate(x,px)
	matplot(px,B%*%y,xlab=xlab,ylab=ylab,type=type,...)
})


OuterProdSecondDerivative<-function(basis){
	M<-basis@Matrices
	k<-basis@order
	d<-dim(M)
	Delta<-1/Hankel(1:(2*k-1),k,k)
	D2<-MatrixPower(DerivativeMatrix(k),2)
	OPSD<-0
	di<-diff(basis@knots[(k):(length(basis@knots)-k+1)])
	for(i in 1:d[3])OPSD<-OPSD+di[i]*t(M[,,i])%*%D2%*%Delta%*%t(D2)%*%M[,,i]
	OPSD
}
fitLS<-function(object, x, y, penalty=0){
	stopifnot(is(object,"SplineBasis"))
	stopifnot(is(x,"numeric"))
	stopifnot(is(y,"numeric"))
	B<-evaluate(object,x)
	solve(crossprod(B)+penalty*OuterProdSecondDerivative(object),crossprod(B,y))
}

setGeneric("penaltyMatrix",function(object, ...)standardGeneric("penaltyMatrix"))
setMethod("penaltyMatrix",signature("SplineBasis"),function(object,...)OuterProdSecondDerivative(object))



#Helper Functions
DerivativeMatrix<-function(n){
	A<-rbind(0,diag(x=1,nrow=n-1,ncol=n))
	B<-diag(1:n,nrow=n,ncol=n)
	A%*%B
}
MatrixPower<-function(A,n){
	n<-as.integer(n)
	stopifnot(n>=1)
	stopifnot(ncol(A)==nrow(A))
	B<-diag(1,nrow(A))
	for(i in seq_len(n))B<-B%*%A
	return(B)
}
I_Sel<-function(i,k,n){
	cbind(
		matrix(0,nrow=k,ncol=i-(k)),
		diag(nrow=k),
		matrix(0,nrow=k,ncol=n-k-i)
		)
}
FindSplineMatrices<-function(knots,k){
n<-length(knots)-1 # n+1= number of knots

D0<-function(i,j,k){
	a<-knots[i]-knots[j]
	b<-knots[j+k-1]-knots[j]
	rtn<-a/b
	rtn[a==0]<-0
	rtn
}
D1<-function(i,j,k){
	a<-knots[i+1]-knots[i]
	b<-knots[j+k-1]-knots[j]
	rtn<-a/b
	rtn[a==0]<-0
	rtn
} 
M<-function(k,i) { 
	if(k==1) return(1);
	if(i<k|i>length(knots)-k) return(matrix(0,k,k));
	rbind(M(k-1,i),0)%*%(cbind(diag(x=1-D0(i,(i-k+2):i,k),k-1),0)+cbind(0,diag(x=D0(i,(i-k+2):i,k),k-1)))+
	rbind(0,M(k-1,i))%*%(cbind(diag(x= -D1(i,(i-k+2):i,k),k-1),0)+cbind(0,diag(x=D1(i,(i-k+2):i,k),k-1)))
} 
MFinal<-array(dim=c(k,k,length(knots)))
for(i in 1:length(knots)){
	MFinal[,,i]<-M(k,i)
}
class(MFinal)<-'SplineMatrices'

MFinal
}
Hankel<-function(x,nrow=length(x)%/%2,ncol=length(x)%/%2){
	Z<-matrix(x[1:ncol],nrow=nrow,ncol=ncol,byrow=T)
	for(i in 1:(nrow-1)){
		Z[i+1,]<-x[-(1:i)][1:ncol]
	}
	Z
}
expand.knots<-function(interior, order=4){
	knots<-interior[c(
	rep(1,order-1),
	seq(length(interior)),
	rep(length(interior),order-1))]
	attr(knots,'order')<-order
	knots
}
#Aliasing
OBasis<-function(...)OrthogonalSplineBasis(...)
