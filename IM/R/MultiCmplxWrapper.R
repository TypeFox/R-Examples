
# Wrapper functions for Multiple Images 
# GPZM, FM, CHF, RF
# Author: Tan Dang
###############################################################################


GPZMmultiple <- function(images=NULL,radius=NULL,theta=NULL,order=NULL,a=NULL, storePoly=FALSE) {
	# double *I, double *ReM, double *ImM, int *dim, double *radius, double *theta, 
	#	int *order, int *a2, double *c, int*NF
	
	d = dim(images[[1]]);
	
	#check input arguments
	if(is.null(images) & is.null(radius) & is.null(theta) & is.null(order) & is.null(a)) {
		sprintf("input arguments: (I=image,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)\n")
		return();
	}
	if(is.null(images)) {
		sprintf("Missing input list of image");
		return();
	}else if(length(d)!=2) {
		sprintf("First argument must be a 2 dimensional image matrix")
		return();
	}
	if(is.null(order)) {
		sprintf("missing maximum order argument - default to 20")
		order=20;
	}
	if(is.null(a)) {
		sprintf("missing stabilization parameter a - default to 0")
		a=0;
	}
	
	nF= length(images);

	# Create M to store moment
	ReM=array(0,c(nF,order+1,order+1));
	ImM=array(0,c(nF,order+1,order+1));
	
	IM= array(0,c(nF,d[1],d[2]));
		
	for (f in 1:nF){
		IM[f,,]=images[[f]];
	}
	
	const= array(0,c(order+1, order+1));
	
	if (storePoly){
		# 4D matrix to store polynomial values - limit to order 12
		polyVals= array(0, c(d[1], d[2], order+1, order+1));
	} else {
		polyVals= array(1);
	}
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.load("src/CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.load("src/CmplxMultiple.dll")
	#}
	
	
	# Compute moment
	out=.C("GPZMMultiplebyC",I=as.double(IM),ReM=as.double(ReM),ImM=as.double(ImM),
			D=as.integer(c(d[1],d[2])),radius=as.double(radius),
			theta=as.double(theta),order=as.integer(order),a=as.integer(a), const=as.double(const),
			NF=as.integer(nF), polyValsR=as.double(polyVals), polyValsI=as.double(polyVals), storePoly=as.integer(storePoly), PACKAGE="IM");
	
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.unload("CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.unload("CmplxMultiple.dll")
	#}
	
	M=complex(length(out$ReM),real=out$ReM,imaginary=out$ImM);
	M=matrix(M,nrow=nF,byrow=FALSE);
	
	moments=list();
	for (f in 1:nF){
		moments= c(moments, list(matrix(M[f,],nrow=(order+1),ncol=(order+1),byrow=FALSE)) );
	}
	
	if (storePoly){
		# Get polyvals
		polyVals= complex(length(out$polyValsR), real=out$polyValsR, imaginary=out$polyValsI);
		polyVals= array(polyVals, c(d[1], d[2], order+1, order+1));
	} else {
		polyVals= NULL;
	}
	list(moments, polyVals)
}


Cmplxmultiple <- function(type,images=NULL,radius=NULL,theta=NULL,nmax=NULL,mmax=NULL, storePoly=FALSE) {
	
	if (type=="fc"){
		run="CHFMMultiplebyC";	
	} else if (type=="fm"){
		run="FMMultiplebyC";
	} else if (type=="fr"){
		run="RFMultiplebyC";
	} else {
		print("Type non-support");
		return(NULL)
	}
	# (double *I, double *ReM, double *ImM,  double *rmat, double *thetamat, 
	# int *D, int *Nmax, int *Mmax, int *NF)
	
	d = dim(images[[1]]);
	
	#check input arguments
	if(is.null(images) & is.null(radius) & is.null(theta) & is.null(nmax) & is.null(mmax)) {
		sprintf("input arguments: (images=images,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)\n")
		return();
	}
	if(is.null(images)) {
		sprintf("Missing input list of image");
		return();
	}else if(length(d)!=2) {
		sprintf("First argument must be a 2 dimensional image matrix")
		return();
	}
	if(is.null(nmax)) {
		sprintf("missing maximum order argument - default to 20")
		nmax=20;
	}
	if(is.null(mmax)) {
		sprintf("missing maximum order argument - default to 20")
		mmax=20;
	}
	
	nF= length(images);
	
	# Create M to store moment
	ReM=array(0,c(nF,nmax+1,mmax+1));
	ImM=array(0,c(nF,nmax+1,mmax+1));
	
	IM= array(0,c(nF,d[1],d[2]));
		
	for (f in 1:nF){
		IM[f,,]=images[[f]];
	}
	
	if (storePoly){
		# store polynomial values - limit to order 12
		polyVals= array(0, c(d[1], d[2], nmax+1, mmax+1));
	} else {
		polyVals = array(1)
	}
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.load("src/CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.load("src/CmplxMultiple.dll")
	#}
	
	# Compute moment
	out=.C(run,I=as.double(IM),ReM=as.double(ReM),ImM=as.double(ImM),radius=as.double(radius),
			theta=as.double(theta),D=as.integer(d),Nmax=as.integer(nmax),Mmax=as.integer(mmax),NF=as.integer(nF), 
			polyValsR=as.double(polyVals), polyValsI=as.double(polyVals), storePoly=as.integer(storePoly), PACKAGE="IM");
		
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.unload("src/CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.unload("src/CmplxMultiple.dll")
	#}
	
	M=complex(length(out$ReM),real=out$ReM,imaginary=out$ImM);
	M=matrix(M,nrow=nF,byrow=FALSE);
	
	moments=list();
	for (f in 1:nF){
		moments= c(moments, list(matrix(M[f,],nrow=(nmax+1),ncol=(mmax+1),byrow=FALSE)) );
	}
	
	if (storePoly){
		polyVals= complex(length.out=length(out$polyValsR),real=out$polyValsR,imaginary=out$polyValsI)
		polyVals= array(polyVals, c(d[1], d[2], nmax+1, mmax+1));
	} else {
		polyVals= NULL;
	}
	list(moments, polyVals)
}


GPZMreconMulti <- function(moments=NULL,radius=NULL,theta=NULL,order=NULL,a=NULL) {
	# double *I, double *ReM, double *ImM, int *dim, double *radius, double *theta, 
	#	int *order, int *a2, double *c, int*NF
	
	d = dim(radius);
	
	#check input arguments
	if(is.null(moments) & is.null(radius) & is.null(theta) & is.null(order) & is.null(a)) {
		sprintf("input arguments: (moments,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)\n")
		return();
	}
	if(is.null(moments)) {
		sprintf("Missing input list of moments");
		return();
	}else if(length(d)!=2) {
		sprintf("First argument must be a 2 dimensional image matrix")
		return();
	}
	if(is.null(order)) {
		sprintf("missing maximum order argument - default to 20")
		order=20;
	}
	if(is.null(a)) {
		sprintf("missing stabilization parameter a - default to 0")
		a=0;
	}
	
	nF= length(moments);

	# Stack Real and Imaginary moments in 1 matrix
	ReM=array(0,c(nF,order+1,order+1));
	ImM=array(0,c(nF,order+1,order+1));

	for (f in 1:nF){
		ReM[f,,]= Re(moments[[f]])[1:(order+1), 1:(order+1)];
		ImM[f,,]= Im(moments[[f]])[1:(order+1), 1:(order+1)];
	}
	
	const= array(0, c(order+1, order+1));
	
	# Empty matrix of images
	IM= array(0,c(nF,d[1],d[2]));
	
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.load("src/CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.load("src/CmplxMultiple.dll")
	#}
	
	
	# Compute moment
	out=.C("GPZMreconMulti",IR=as.double(IM),ReM=as.double(ReM),ImM=as.double(ImM),D=as.integer(c(d[1],d[2])),
			radius=as.double(radius), theta=as.double(theta),Pmax=as.integer(order),a1=as.integer(a),
			const=as.double(const), NF=as.integer(nF), PACKAGE="IM");
	
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.unload("CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.unload("CmplxMultiple.dll")
	#}
	
	
	IM=matrix(out$IR,nrow=nF,byrow=FALSE);

	images=list();
	for (f in 1:nF){
		#images= c(images, list(matrix(IM[f,],nrow=(d[1]),ncol=(d[2]),byrow=FALSE)) );
		images= c(images, list(flip(matrix(IM[f,],nrow=(d[1]),ncol=(d[2]),byrow=FALSE))) );
	}
	
	images
}


CmplxreconMulti <- function(type, moments=NULL,radius=NULL,theta=NULL,nmax=NULL,mmax=NULL) {
	if (type=="fc"){
		run="CHFreconMulti";	
	} else if (type=="fm"){
		run="FMreconMulti";
	} else if (type=="fr"){
		run="RFreconMulti";
	} else {
		print("Type non-support");
		return(NULL)
	}
	
	# double *I, double *ReM, double *ImM, int *dim, double *radius, double *theta, 
	#	int *order, int *a2, double *c, int*NF
	
	d = dim(radius);
	
	#check input arguments
	if(is.null(moments) & is.null(radius) & is.null(theta) & is.null(nmax) & is.null(mmax)) {
		sprintf("input arguments: (moments,radius=radius of polar representation, 
						theta=angle of polar representation, order=maximum moment order,
						a=stabilization parameter)\n")
		return();
	}
	if(is.null(moments)) {
		sprintf("Missing input list of moments");
		return();
	}else if(length(d)!=2) {
		sprintf("First argument must be a 2 dimensional image matrix")
		return();
	}
	if(is.null(nmax)) {
		sprintf("missing maximum order argument - default to 20")
		nmax=20;
	}
	if(is.null(mmax)) {
		sprintf("missing maximum order argument - default to 20")
		mmax=20;
	}
	
	nF= length(moments);

	# Stack Real and Imaginary moments in 1 matrix
	ReM=array(0,c(nF,nmax+1,mmax+1));
	ImM=array(0,c(nF,nmax+1,mmax+1));

	for (f in 1:nF){
		ReM[f,,]= Re(moments[[f]])[1:(nmax+1), 1:(mmax+1)];
		ImM[f,,]= Im(moments[[f]])[1:(nmax+1), 1:(mmax+1)];
	}
	
	# Empty matrix of images
	IM= array(0,c(nF,d[1],d[2]));
	
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.load("src/CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.load("src/CmplxMultiple.dll")
	#}
	
	
	# Compute moment
	out=.C(run,IR=as.double(IM),ReM=as.double(ReM),ImM=as.double(ImM),radius=as.double(radius),
			theta=as.double(theta),D=as.integer(d),Nmax=as.integer(nmax),Mmax=as.integer(mmax),NF=as.integer(nF), PACKAGE="IM");
	
	#check operating system
	#if(.Platform$OS.type == "unix") {
	#	dyn.unload("CmplxMultiple.so")
	#}else if(.Platform$OS.type == "windows") {
	#	dyn.unload("CmplxMultiple.dll")
	#}
	
	
	IM=matrix(out$IR,nrow=nF,byrow=FALSE);

	images=list();
	for (f in 1:nF){
		#images= c(images, list(matrix(IM[f,],nrow=(d[1]),ncol=(d[2]),byrow=FALSE)) );
		images= c(images, list(flip(matrix(IM[f,],nrow=(d[1]),ncol=(d[2]),byrow=FALSE))));
	}
	
	images
}