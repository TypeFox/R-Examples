# wrapper functions for gpzm.cpp 
# 
# Author: Allison Irvine, Tan Dang
###############################################################################


GPZMmoments <- function(I=NULL,radius=NULL,theta=NULL,order=NULL,a=NULL) {
	
	errors=""
	M=0
	d = dim(I);
	
	#check input arguments
	if(is.null(I) & is.null(radius) & is.null(theta) & is.null(order) & is.null(a)) {
		errors = "input arguments: (I=image,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)"
		return(list(M,errors));
	}
	if(is.null(I)) {
		errors = "Missing input image";
		return(list(M,errors));
	}else if(length(d)!=2) {
		errors = "First argument must be a 2 dimensional image matrix"
		return(list(M,errors));
	}
	if(is.null(order)) {
		errors = "Missing maximum order argument"
		return(list(M,errors));
	}
	if(is.null(a)) {
		errors = "Missing stabilization parameter a - default to 0"
		a=0;
	}
	if(length(order)>1) {
		errors = "repetition is constrained by order - using first value in the variable order"
		order = order[1];
	}
	
	#check operating system
	#if(.Platform == "unix") {
	#	dyn.load("cppFiles/gpzm.so")
	#}else if(.Platform == "windows") {
	#	dyn.load("cppFiles/gpzm.dll")
	#}
	
	#gpzmMoments(double *I, double *ReM, double *ImM ,int *dim, double *radius, double *theta, int *order,
	#double *a2, double *c, double *polyValsR, double *polyValsI, int *storePoly)
	
	#call C function
	out=.C("gpzmMoments",I=as.double(I),ReM=as.double(mat.or.vec(order+1,order+1)),
			ImM=as.double(mat.or.vec(order+1,order+1)),
			d=as.integer(d),radius=as.double(radius),theta=as.double(theta), 
			order=as.integer(order), a=as.double(a), c=as.double(mat.or.vec(order+1,order+1)), PACKAGE="IM");
	
	M = complex(length.out=length(out$ReM),real=out$ReM,imaginary=out$ImM)
	M = matrix(M,nrow=order+1,ncol=order+1,byrow=FALSE)
	constants = matrix(out$c,nrow=order+1,ncol=order+1,byrow=FALSE)
	
	list(M,constants,errors)
	
}

#reconstruct image
GPZMrecon <- function(M=NULL,radius=NULL,theta=NULL,order=NULL,a=NULL,constants=NULL) {
	
	errors = ""
	IR=0
	d = dim(radius)
	
	#check input arguments
	if(is.null(M) & is.null(radius) & is.null(theta) & is.null(order) & is.null(a)) {
		errors = "input arguments: (M=image moments,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)"
		return(list(IR,errors));
	}
	if(is.null(radius)) {
		errors = "Missing pixel radius"
		return(list(IR,errors));
	}	
	if(is.null(M)) {
		errors = "Missing image moments"
		return(list(IR,errors));
	}else if(length(d)!=2) {
		errors = "First argument must be a 2 dimensional matrix"
		return(list(IR,errors));
	}
	if(is.null(order)) {
		errors = "Missing maximum order argument - default to maximum available"
		order = checkOrder(M);
		order = order[1];
	}
	if(is.null(a)) {
		errors = "Missing stabilization parameter a - default to 0"
		a=0;
	}
	
	#check operating system
	#if(.Platform == "unix") {
	#	dyn.load("cppFiles/gpzm.so")
	#}else if(.Platform == "windows") {
	#	dyn.load("cppFiles/gpzm.dll")
	#}
	
	#gpzmReconstruct(double *IR, complex <double> *M, int *d, double *radius, double *theta, 
	#		int *order, double *a2, double *c) 

	#call C function
	out = .C("gpzmReconstruct",IR=as.double(mat.or.vec(d[1],d[2])), ReM=as.double(Re(M[1:(order[1]+1),1:(order[1]+1)])), ImM=as.double(Im(M[1:(order[1]+1),1:(order[1]+1)])),d=as.integer(dim(radius)), 
			radius=as.double(radius),theta=as.double(theta), order=as.integer(order),
			a2=as.double(a), c=as.double(constants), PACKAGE="IM")
	
	IR = matrix(out$IR,nrow=d[1],ncol=d[2],byrow=FALSE)
	
	return(list(flip(IR),errors));
}
