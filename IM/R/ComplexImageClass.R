#ComplxIm Class
# inherits from Image class, contains additional parameters required for complex moments analysis
# 
# Author: Allison Irvine, Tan Dang
###############################################################################

#define CmplxIm class for complex moments calculations
setClass(
	Class="CmplxIm",
	representation = representation(
			radius = "matrix",
			theta = "matrix",
			constants = "matrix",
			params = "numeric",
			order = "numeric"
	),
	contains="Image"
)

#define initialization with an input image or filename
#define initialization
setMethod(
		f="initialize",
		signature="CmplxIm",
		definition = function(.Object,img=NULL,filename = "") {
			#initialize inherited attributes
			if((!is.null(img)) && ((dim(img)[1]>0) && (dim(img)[2]>0))) {
				as(.Object,"Image") <- new("Image",img=img);
			} else if(filename != "") {
				as(.Object,"Image") <- new("Image",filename=filename);
			}else{
				.Object@error = "No input image"
				return(.Object)
			}
			
			#set pixel radius and theta values
			setPolar(.Object) = calcCentroid(.Object@I);
			return(.Object)
			
		}
)

########### SETTER FUNCTIONS ####################
setReplaceMethod(
		f="setImage",
		signature=c("CmplxIm"),
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
			if((dim(obj@I)[1]>0) && (dim(obj@I)[2]>0)){
				#set pixel radius and theta values
				setPolar(obj) = NULL;
			}
			return(obj)
		}
)

#function definition to set polar coordinates of pixels - can pass a centroid value, if null, it is calculated
setGeneric (
		name = "setPolar<-",
		def=function(obj,value) {standardGeneric("setPolar<-")},
)

setReplaceMethod(
		f="setPolar",
		signature=c("CmplxIm"),
		definition = function(obj,value=NULL) {
			if(is.null(value)) {
				#calculate centroid
				obj@centroid = calcCentroid(obj@I);
			}else if(length(value) != 2) {
				obj@error[length(obj@error)+1] = "centroid must be a numeric array of length 2 - automatically calculating centroid"
				obj@centroid = calcCentroid(obj@I);
			}else{
				obj@centroid = value
			}
			
			p = polarXY(obj@I,obj@centroid);
			obj@radius = p[[1]];
			obj@theta = p[[2]];
			return(obj)
		}
)

#set moment type
#choice of:
#gpzm		Generalized Pseudo-Zernike
#fm			fourier mellin
#fc			fourier chebyshev
#fr			radial fourier
setReplaceMethod(
		f="momentType",
		signature=c("CmplxIm"),
		definition = function(obj,value) {
			#make sure it is one of the available
			if( (sum(value == c("gpzm","fm","fr","fc")) != 1) || (length(value)!=1) ){
				obj@error[length(obj@error)+1] = "Please choose one of the following moment types:gpzm,fm,fr,fc"
				return(obj)
			}
			
			obj@momentType = value
			return(obj)
		}
)




setReplaceMethod(
		f="setOrder",
		signature=c("CmplxIm"),
		definition = function(obj,value) {
			
			#check input values
			if(!is.numeric(value)) {
				obj@error[length(obj@error)+1] = "order must be a numeric array of length 2"
				return(obj)
			}

			#set params for individual moment calculations
			if(obj@momentType=="gpzm") {	#generalized pseudo-zernike
				if(length(value)>1) {
					obj@error[length(obj@error)+1] = "Generalized pseudo-zernike moments only require one order value. Repetition is constrained by order. Only first value being used."
					obj@order = value[1]
					return(obj)
				}
				obj@order = value
			} else if (length(value)!=2){
				obj@error[length(obj@error)+1] = "order must be a numeric array of length 2: c(order,repetition)"
				return(obj)
			}else {
				obj@order = value
			}
			
			return(obj)
		}
)

setReplaceMethod(
		f="setParams",
		signature=c("CmplxIm"),
		definition = function(obj,value) {
			if(obj@momentType!="gpzm") {
				obj@error[length(obj@error)+1] = "only generalized pseudo-zernike moments require a parameter"
				return(obj)
			}
			if(!is.numeric(value) || length(value)!=1) {
				obj@error[length(obj@error)+1] = "only one parameter is required for generalized pseudo-zernike moments"
				return(obj)
			}
			
			obj@params = value
			return(obj)		
		}
)


################# AFTER INITIALIZATION ##########################

#function definition for moments calculation
setReplaceMethod(
		f="Moments",
		signature=c("CmplxIm"),
		definition = function(obj,value) {
			
			#call wrapper functions
			if(obj@momentType=="gpzm") {			#generalized pseudo-zernike
				result = GPZMmoments(obj@I,obj@radius,obj@theta,obj@order,obj@params)
				obj@moments = result[[1]];
				obj@constants = result[[2]];
				if (result[[3]] != "") {
					obj@error[length(obj@error)+1] = result[[3]];
				}
				return(obj)
			} else if(obj@momentType == "fm") {		#fourier mellin
				result = FMmoments(obj@I,obj@radius,obj@theta,obj@order)
			} else if(obj@momentType=="fr") {		#radial fourier
				result = FRmoments(obj@I,obj@radius,obj@theta,obj@order)
			} else if(obj@momentType=="fc") {		#fourier chebyshev
				result = FCmoments(obj@I,obj@radius,obj@theta,obj@order)
			} else {
				obj@error[length(obj@error)+1] = "Moment type not specified."
				return(obj)
			}
			
			obj@moments = result[[1]];
			if (result[[2]] != "") {
				obj@error[length(obj@error)+1] = result[[2]];
			}
			
			return(obj)
		}
)



#function definition for image reconstruction from moments
#value is c(order,repetition) 
setReplaceMethod(
		f="Reconstruct",
		signature=c("CmplxIm"),
		definition = function(obj,value) {
			
			#if no maximum order given for reconstruction, use maximum available (not nan or inf)
			order = checkOrder(obj@moments);

			#if a maximum order is provided, check its validity
			if(!is.null(value)) {
				if(value[1]>order[1]) {
					obj@error[length(obj@error)+1] = sprintf("maximum order is %d. Order is being reset to this value.",order[1])
					value[1] = order[1]
				}
				if(obj@momentType != "gpzm") {
					if(value[2]>order[2]) {
						obj@error[length(obj@error)+1] = sprintf("maximum repetition is %d. Repetition is being reset to this value.",order[2])
						value[2] =  order[2]
					}
				}
				order=value
			}else{
				obj@error[length(obj@error)+1] = "No order and repetition specified. Using maximum available."
			}
			
			#call wrapper functions
			if(obj@momentType=="gpzm") {			#generalized pseudo-zernike
				result = GPZMrecon(obj@moments,obj@radius,obj@theta,min(order),
						obj@params,obj@constants)
			} else if(obj@momentType == "fm") {		#fourier mellin
				result = FMrecon(obj@moments,obj@radius,obj@theta,order)
			} else if(obj@momentType=="fr") {		#radial fourier
				result = FRrecon(obj@moments,obj@radius,obj@theta,order)
			} else if(obj@momentType=="fc") {		#fourier chebyshev
				result = FCrecon(obj@moments,obj@radius,obj@theta,order)
			}else{
				obj@error[length(obj@error)+1] = "Please choose one of the following moment types:\ngpzm,fm,fr,fc"
				return(obj)
			}
				
			obj@reconstruction = result[[1]];
			if(result[[2]]!=""){
				obj@error[length(obj@error)+1] = result[[2]];
			}
			return(obj)
		}
)

setReplaceMethod(
		f="Invariant",
		signature= c("CmplxIm"),
		definition = function(obj, value){
			#if moments have not already been calculated, do so now
			if (sum(dim(obj@moments))==0){
				Moments(obj)= NULL;
			}
			M= abs(obj@moments);
			obj@invariant=M/M[1,1];
			return(obj)
		}
)