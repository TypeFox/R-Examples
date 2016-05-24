# OrthIm Class
# inherits from Image class, contains additional parameters required for discrete/continuous orthogonal moments analysis
# 
# Author: Allison Irvine, Tan Dang
###############################################################################

#define OrthIm class for discrete/continuous real orthogonal moments calculations
setClass(
		Class="OrthIm",
		representation = representation(
				polynomials = "list",
				params = "list",
				order = "numeric"
		),
		contains="Image"
)

#define initialization with input image or filename
setMethod(
		f="initialize",
		signature="OrthIm",
		definition = function(.Object,img=NULL,filename="") {
			
			#initialize inherited attributes
			if((!is.null(img)) && ((dim(img)[1]>0) && (dim(img)[2]>0))) {
				as(.Object,"Image") <- new("Image",img=img);
			} else if(filename != "") {
				as(.Object,"Image") <- new("Image",filename=filename);
			}else{
				.Object@error[length(obj@error)+1] = "No input image"
				return(.Object)
			}

			return(.Object)

		}
)



################# AFTER INITIALIZATION ##########################
#function definition to set moments type and parameters

####### can choose a combination of 2 from the same group: 
#DISCRETE:
#"cheby"		Chebyshev
#"krawt"		Krawtchouk
#"hahn"			Hahn

#CONTINUOUS:
#"gegen"		Gegenbauer
#"chebycont"	ChebyshevCont	(continuous chebyshev)
#"legend"		Legendre

setReplaceMethod(
		f="momentType",
		signature=c("OrthIm"),
		definition = function(obj,value) {
			#check that up to 2 types of polynomials were chosen
			if(length(value) > 2) {
				obj@error[length(obj@error)+1] = "You can only choose up to 2 different polynomials to combine"
				return(obj)
			}else if(length(value)>1) {
				#check that we are combining 2 of a matching type (both continuous or discrete)
				#if(( ( (sum(value[1]==c("cheby","krawt","hahn"))
				#				+ (sum(value[2]==c("cheby","krawt","hahn")) )%%2)==1)))  {
				#	sprintf("The pair of polynomial types must both be either discrete or continuous")
				#	return(obj)
				#}
			}else if(length(value)==1) {
				value[2] = value[1]
			}
			
			#check that the moment type is one of the available
			if((sum(value[1] == c("cheby","krawt","hahn","gegen","chebycont","legend"))
						+sum(value[2] == c("cheby","krawt","hahn","gegen","chebycont","legend")))<2) {
				obj@error[length(obj@error)+1] = "you must choose 1 or 2 of the available moment types:"
				return(obj)
			}
			
			obj@momentType = value;
			return(obj)
		}
)



#value is a list of length 2
#the elements of the list are each a numeric array of parameters, one each for the x and y directions
#params[[1]] = c(a,b)
setReplaceMethod(
		f="setParams",
		signature=c("OrthIm"),
		definition = function(obj,value) {
			
			tempValue = value;
			
			#check input values
			if(!is.null(dim(value)) || length(value) > 2) {
				obj@error[length(obj@error)+1] = "params must be a list of length 2: parameters for x, parameters for y"
				return(obj)
			}
			
			if(is.null(value)) {
				obj@error[length(obj@error)+1] = "parameters not set"
				if(length(obj@params)==0) {
					obj@params = list(NULL,NULL);
				}
			}
			
			if (length(value)==1){
				tempValue = c(as.list(value), as.list(value));	
			} 
			
			#value=as.list(value)
			#set param names, more error checking
			for (i in 1:2) {
				type = obj@momentType[i]; 
				
				
				if(type=="cheby" || type=="chebycont" || type=="legend") {	#Chebyshev
					obj@params[[i]] = NULL;
				} else if((type == "krawt")||(type=="gegen")) {		#Krawtchouk,Gegenbauer
					if(length(value[[i]]) < 1) {
						obj@error[length(obj@error)+1] = sprintf("%s missing parameters: 1 parameter required: a",type)
						tempValue[[i]]=NULL;
					} else if(length(value[[i]]) == 1) {
						#check parameter constraints for krawtchouk
						if ( (type=="krawt") && ((value[[i]]<0) || (value[[i]]>1)) ){
							obj@error[length(obj@error)+1] = "constraint (0<a<1) is not satisfied"
							tempValue[[i]]=NULL;
							return(obj)
						}
						#check parameter constraints for gegenbauer
						if ( (type=="gegen") && ((value[[i]]<(-1/2)) || (value[[i]]==0)) ){
							obj@error[length(obj@error)+1] = "constraints (a != 0) and (a > -1/2) are not satisfied"
							tempValue[[i]]=NULL;
							return(obj)
						}
						names(tempValue[[i]]) = c("a")
					} else {
						obj@error[length(obj@error)+1] = sprintf("%s parameters: 1 parameter required: a",type)	
						tempValue[[i]]=NULL;
					} 
				} else if (type=="hahn") {	#Hahn
					if((length(value[[i]]) < 2) || is.null(value[[i]])) {
						obj@error[length(obj@error)+1] = sprintf("%s missing parameters: 2 parameters required: a, c",type)
						tempValue[[i]]=NULL;
					} else if (length(value[[i]]) == 2) {
						if (value[[i]][1] <= (-1/2)) {
							obj@error[length(obj@error)+1] = "constraint (a > -1/2) is not satisfied"
							tempValue[[i]]=NULL;
							return(obj)	
						}
						if (value[[i]][1] <= (abs(value[[i]][2])-1)){
							obj@error[length(obj@error)+1] = "constraint (a> abs(c)-1) is not satisfied"
							tempValue[[i]]=NULL;
							return(obj)
						} 
						#if all conditions are met,
						names(tempValue[[i]]) = c("a","c") 
					} else {
						obj@error[length(obj@error)+1] = sprintf("%s: 2 parameters required: a, c",type)
						tempValue[[i]]=NULL;
					} 
				}else{
					obj@error[length(obj@error)+1] = "Set the moment type first > momentType(imageObject) = momentType"
					return(obj)
				}
			}
			
			#bug fix
			if (length(tempValue)==1) {
				tempValue = list(tempValue[[1]],NULL);
			}
			
			#set params for individual moment calculations
			obj@params = tempValue;
			names(obj@params)=c("Xparameters","Yparameters")
			
			return(obj)
		}
)
			




#value is a numeric array of length 2
#the elements of the list are the orders for the x and y directions
setReplaceMethod(
		f="setOrder",
		signature=c("OrthIm"),
		definition = function(obj,value) {
			#check input values
			if(!is.null(dim(value)) || (length(value) > 2) || (length(value) < 1)) {
				obj@error[length(obj@error)+1] = "params must be a numeric array of length 2: order for x, order for y"
				return(obj)
			}
			if (length(value)==1){
				value= c(value, value);	
			}
			obj@order = value;
			return(obj)
			
		}
)
			


#function definition for moments calculation
setReplaceMethod(
		f="Moments",
		signature=c("OrthIm"),
		definition = function(obj,value) {
			#check for neccessary parameters
			if (( (sum(obj@momentType[1]==c("hahn","krawt","gegen"))==1) + sum(obj@momentType[2]==c("hahn","krawt","gegen")) ) > length(obj@params)) {
				obj@error[length(obj@error)+1] = "Missing parameters"
			} else if ((sum(obj@momentType[1]==c("hahn","krawt","gegen"))==1) && is.null(obj@params[[1]])) {
				obj@error[length(obj@error)+1] = sprintf("Missing parameters for %s  polynomials",obj@momentType[1]);
			} else if ((sum(obj@momentType[2]==c("hahn","krawt","gegen"))==1) && is.null(obj@params[[2]])) {
				obj@error[length(obj@error)+1] = sprintf("Missing parameters for %s  polynomials",obj@momentType[2]);
			} else {
				#calculate polynomials
				obj@polynomials= OrthPolynomials(obj@momentType,obj@dimensions,obj@order,obj@params);
				#calculate moments
				obj@moments = OrthMoments(obj@I,obj@polynomials);
			}
			
			return(obj)
		}
)

#function definition for moments calculation
setReplaceMethod(
		f="Invariant",
		signature=c("OrthIm"),
		definition = function(obj,value) {
			storeObj= obj;
			
			I=polarTransform(obj@I,value[[1]],value[[2]], calcCentroid(obj@I))[[1]];
			
			setImage(obj)=I
			
			#calculate moments
			Moments(obj)=NULL;
			
			storeObj@invariant = obj@moments;
			return(storeObj)
		}
)

#function definition for image reconstruction from moments
setReplaceMethod(
		f="Reconstruct",
		signature=c("OrthIm"),
		definition = function(obj,value=NULL) {
			#if no maximum order given for reconstruction, use maximum available (not nan or inf)
			order = checkOrder(obj@moments);
			#if order defined, check values
			if(!is.null(value)) {
				#check that user-defined order is valid
				if(value[1]>order[1]) {
					obj@error[length(obj@error)+1] = sprintf("maximum order in x direction is %d. Order is being reset to this value.",order[1])
					value[1] = order[1]
				}
				if(value[2]>order[2]) {
					obj@error[length(obj@error)+1] = sprintf("maximum order in y direction is %d. Order is being reset to this value.",order[2])
					value[2] =  order[2]
				}
				order=value
			}

			obj@reconstruction = OrthReconstruct(obj@moments,obj@polynomials,obj@dimensions,order)
			
			return(obj)
		}
)


setGeneric (
		name = "transform<-",
		def=function(obj,value) {standardGeneric("transform<-")},
)

#transform image to polar plot of radius vs theta
setReplaceMethod(
		f="transform",
		signature=c("OrthIm"),
		definition = function(obj,value=NULL) {
			resolution = value
			#check input values
			if(is.null(value)) {
				resolution = 8;
				obj@error[length(obj@error)+1] = "no resolution given, defaulting to 8"
			}
			
			I=polarTransform(obj@I,resolution)[[1]];
			
			setImage(obj)= I;
			return(obj)
			
		}
)


setMethod(
		f="plotPoly",
		signature=c("OrthIm","numeric"),
		definition = function(obj,order) {
			P=obj@polynomials;
			
			#plot up to 20 different colors
			colorVec = c(10,12,17,24,26,31,32,33,43,47,51,54,60,62,73,75,79,91,142,259) 
			typeName = ""

			nplots=1;
			if(obj@momentType[1]!=obj@momentType[2]) {
				par(mfrow = c(2,1))
				nplots=2;
			}
					
			for (i in 1:nplots) {
				#get moment type full name
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
				
				#if a specific order(s) was given, plot  that order only
				if (!is.null(order)) {
					X = P[[i]][,order[1]];
					plot(1:length(X),X,type="n",xlab="X coordinate",ylab="polynomial(X)")
					title(main = sprintf("%s Polynomials, order %d to %d",typeName[i],order[1],order[length(order)]))
					h = 1;
					for(j in order) {
						X = P[[i]][,j];
						lines(1:length(X),X,col=colors()[colorVec[h%%20]])
						h=h+1
					}
				}
				
			}
			
		}
)
