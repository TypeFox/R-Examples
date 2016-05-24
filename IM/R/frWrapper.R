# wrapper functions for RadialFourier.cpp 
# 
# Author: Allison Irvine, Tan Dang
###############################################################################


FRmoments <- function(I=NULL,radius=NULL,theta=NULL,order=NULL) {
	
	errors = ""
	M = 0;
	d = dim(I);
	
	#check input arguments
	if(is.null(I) & is.null(radius) & is.null(theta) & is.null(order)) {
		errors = "input arguments: (I=image,radius=radius of polar representation, 
				theta=angle of polar representation, order=maximum moment order)"
		return(list(M,errors));
	}
	if(is.null(I)) {
		errors = "Missing input image"
		return(list(M,errors));
	}else if(length(d)!=2) {
		errors = "First argument must be a 2 dimensional image matrix"
		return(list(M,errors));
	}
	if(is.null(radius)) {
		errors = "Missing radius of polar pixel coordinates"
		return(list(M,errors));
	}
	if(is.null(theta)) {
		errors = "Missing angle of polar pixel coordinates"
		return(list(M,errors));
	}
	if(is.null(order)) {
		errors = "Missing maximum order argument"
		return(list(M,errors));
	}
	
	#check operating system
	#if(.Platform == "unix") {
	#	dyn.load("cppFiles/RadialFourier.so")
	#}else if(.Platform == "windows") {
	#	dyn.load("cppFiles/RadialFourier.dll")
	#}
	
	#RFbyC(double *I, double *ReM, double *ImM,  
	#		double *radius, double *theta, int *D1, int *D2, int *Nmax, int *Mmax){
	
	#call C function
	out=.C("RFbyC",I=as.double(I),ReM=as.double(mat.or.vec(order[1]+1,order[2]+1)),
			ImM=as.double(mat.or.vec(order[1]+1,order[2]+1)),radius=as.double(radius),
			theta=as.double(theta),D=as.integer(d),
			Nmax=as.integer(order[1]), Mmax=as.integer(order[2]), PACKAGE="IM");
	
	M = complex(length.out=length(out$ReM),real=out$ReM,imaginary=out$ImM)
	M = matrix(M,nrow=order[1]+1,ncol=order[2]+1,byrow=FALSE)
	
	list(M,errors)	
}

#reconstruct image
FRrecon <- function(M=NULL,radius=NULL,theta=NULL,order=NULL) {
	
	errors = ""
	IR = 0;
	d=dim(radius);
	
	#check input arguments
	if(is.null(M) & is.null(radius) & is.null(theta) & is.null(order)) {
		errors = "input arguments: (M=image moments,radius=radius of polar representation, 
				theta=angle of polar representation, order=maximum moment order)"
		return(list(IR,errors));
	}
	if(is.null(I)) {
		errors = "Missing moments"
		return(list(IR,errors));
	}else if(length(dim(M))!=2) {
		errors = "First argument must be a 2 dimensional moment matrix"
		return(list(IR,errors));
	}
	if(is.null(radius)) {
		errors = "Missing radius of polar pixel coordinates"
		return(list(IR,errors));
	}
	if(is.null(theta)) {
		errors = "Missing angle of polar pixel coordinates"
		return(list(IR,errors));
	}
	if(is.null(order)) {
		errors = "Missing maximum order argument - default to maximum available"
		order = checkOrder(M)
	}
	
	#check operating system
	#if(.Platform == "unix") {
	#	dyn.load("cppFiles/RadialFourier.so")
	#}else if(.Platform == "windows") {
	#	dyn.load("cppFiles/RadialFourier.dll")
	#}
	
	
	#RFReconstructbyC(double *IR, double *ReM, double *ImM, double *radius, 
	#double *theta, int *D1, int *D2, int *Nmax, int *Mmax){

	#call C function
	out = .C("RFReconstructbyC",IR=as.double(mat.or.vec(d[1],d[2])), ReM=as.double(Re(M[1:(order[1]+1),1:(order[2]+1)])),ImM=as.double(Im(M[1:(order[1]+1),1:(order[2]+1)])),
			radius=as.double(radius),theta=as.double(theta), 
			D=as.integer(d),
			Nmax=as.integer(order[1]), Mmax=as.integer(order[2]), PACKAGE="IM");
	
	IR = matrix(out$IR,nrow=d[1],ncol=d[2],byrow=FALSE)
	
	list(flip(IR),errors)
}
