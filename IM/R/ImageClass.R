#Image Class: S4 class to represent images
# 
# Author: Allison Irvine, Tan Dang
###############################################################################
#requires libraries: jpeg,png,bmp


#basic image class
setClass (
		Class = "Image",
		representation = representation(
				I = "matrix",
				dimensions = "numeric",
				centroid = "numeric",
				filename = "character",
				imType = "character",
				momentType = "character",
				moments = "matrix",
				reconstruction = "matrix",
				invariant = "matrix",
				error = "character"
		)
)

#initialize method, called when a new object is created, requires a filename, if not I@I is empty
setMethod(
		f="initialize",
		signature="Image",
		definition = function(.Object,img=NULL,filename=NULL) {
			
			if(!is.null(img)) {				#create object with input image
				setImage(.Object) <- img;
				if((dim(.Object@I)[1]>0) && (dim(.Object@I)[2]>0)){
					setCentroid(.Object) = NULL;
					.Object@dimensions = dim(.Object@I);
				}else{
					.Object@error = "No input image"
				}
			} else if(!(is.null(filename))) {		#create object with input filename
				.Object@filename=filename;
				.Object@imType="";
				
				#get file type extension
				temp = unlist(strsplit(filename,"\\."));
				
				if(length(temp) > 1) {
					.Object@imType = temp[2];
				}else{
					.Object@error = "filename is not valid";
					return(.Object)
				}
			
				if(.Object@imType=="jpg") {
					I = readJPEG(filename);
				}else if(.Object@imType=="png") {
					I = readPNG(filename);
				}else if(.Object@imType=="bmp") {
					I = read.bmp(filename);
				}else {
					.Object@error = "Not a supported file type"
				}
			
				setImage(.Object) <- img;
			} else {
				.Object@error = "No input image or filename"
			}
			
			return(.Object)
		}		
)


####### SET FUNCTIONS ###########

#set image matrix
setGeneric (
		name= "setImage<-",
		def=function(obj,value) {standardGeneric("setImage<-")},
)

setReplaceMethod(
		f="setImage",
		signature=c("Image"),
		definition = function(obj,value=NULL) {
			
			if(is.null(value)) {
				obj@error[length(obj@error)+1] = "No input image"
				return(obj)
			}else if(length(dim(value))==2) {
				img = value;
			}else if (length(dim(value)) > 2) {
				#convert color image to grayscale
				img = rowSums(value,dims=2);
			}else {
				obj@error[length(obj@error)+1] = "No input image"
				return(obj)
			}
			
			obj@I = as.matrix(img);
			obj@dimensions = dim(obj@I);
			setCentroid(obj) =  NULL;
			return(obj)
		}
)

#calculate centroid
setGeneric (
		name = "setCentroid<-",
		def=function(obj,value) {standardGeneric("setCentroid<-")},
)

setReplaceMethod(
		f="setCentroid",
		signature=c("Image"),
		definition = function(obj,value=NULL) {
			#if no argument passed, centroid is calculated internally
			if(is.null(value)) {
				obj@centroid = calcCentroid(obj@I)
			}else if (length(dim(value))==2){
				obj@centroid = value;
			}else {
				obj@error = "The centroid must be an array of length 2"
			}
			return(obj)
		}
)

#display moments
setGeneric (
		name = "plotMoment",
		def=function(obj) {standardGeneric("plotMoment")},
)

setMethod(
		f="plotMoment",
		signature=c("Image"),
		definition = function(obj) {
			M=obj@moments;
			
			if(is.complex(M)) {
				M=abs(M);
			}
			
			ntypes=length(unique(obj@momentType));
			
			typeName=""
			#get moment type full names
			for(i in 1:ntypes) {
				typeName[i] = switch(obj@momentType[i],
						gpzm = "Generalized Pseudo-Zernike",
						fc = "Fourier Chebyshev",
						fr = "Radial Fourier",
						fm = "Fourier Mellin",
						cheby = "Discrete Chebyshev",
						chebycont = "Continuous Chebyshev",
						legend = "Legendre",
						gegen = "Gegenbauer",
						krawt = "Krawtchouk",
						hahn = "Hahn"
						)
			}
			
			xlab="order in X";
			ylab="order in Y";
			#if moment types are complex
			if (sum(obj@momentType[1]==c("gpzm","fr","fc","fm"))>0) {
				xlab = "order";
				ylab = "repetition";
			}
			
			#plot scaled values
			#rotate image so that it appears aligned
			im <- as.data.frame(t(M));
			im <- rev(im);
			im <- as.matrix(im);
			par(mfrow = c(1,1))
			image((1:dim(im)[1]),(1:dim(im)[2]),asinh(im),xlab=xlab,ylab=ylab)
			if(ntypes==1){
				title(main = sprintf("%s Moments",typeName))
			}else{
				title(main = sprintf("%s - %s Moments",typeName[1],typeName[2]))
			}
			
		}
)


#generic function to set moment type
setGeneric (
		name = "momentType<-",
		def=function(obj,value) {standardGeneric("momentType<-")},
)

#generic function for moments calculation
setGeneric (
		name= "Moments<-",
		def=function(obj,value) {standardGeneric("Moments<-")},
)

#generic function for image reconstruction
setGeneric (
		name= "Reconstruct<-",
		def=function(obj,value) {standardGeneric("Reconstruct<-")},
)

#generic function for computing invariants
setGeneric (
		name= "Invariant<-",
		def=function(obj,value) {standardGeneric("Invariant<-")},
)

#generic function to set order
setGeneric (
		name = "setOrder<-",
		def=function(obj,value) {standardGeneric("setOrder<-")},
)

#generic function to set type-specific moment parameters
setGeneric (
		name = "setParams<-",
		def=function(obj,value) {standardGeneric("setParams<-")},
)

#generic function to display polynomials of orthogonal moments
setGeneric (
		name = "plotPoly",
		def=function(obj,order) {standardGeneric("plotPoly")},
)