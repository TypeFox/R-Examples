# Multiple Images Class
# 
# Author: Allison Irvine, Tan Dang
setClass (
		Class = "MultiIm",
		representation = representation(
				imageList= "list",
				dimension= "numeric",
				polynomials = "list",
				storePoly = "logical",
				momentType = "character",
				order = "numeric",
				Params = "list",
				moments = "list",
				reconstruction = "list",
				invariant = "list",
				error = "character"
		)
)

# initialize method, called when a new object is created
# compute necessary thing, and compute moment	
setMethod(
		f="initialize",
		signature="MultiIm",
		definition = function(.Object,images) {
			
			checkSize(.Object)=images;
			
			.Object@storePoly= FALSE;
			
			return(.Object)
		}
)

# Generics for checkSize method 
setGeneric (
		name= "checkSize<-",
		def=function(obj,value) standardGeneric("checkSize<-"),
)

# Method to check all images are of the same size
setReplaceMethod (
		f="checkSize",
		signature=c("MultiIm"),
		definition = function(obj,value){
			images= value;
			check=TRUE;
			images=lapply(images, function(z){
				if (length(dim(z))>2) {
					return(1/3*(z[,,1]+z[,,2]+z[,,3]))
				} else {
					return(z)
				}
			})
			d=dim(images[[1]]);
			noUse=lapply(images,function(z){
				if (sum(dim(z)!=d)>0) check <<-FALSE; # The <<- assign is needed
			});
			
			if (!check) {
				obj@error[length(obj@error)+1] = "Matrices are not of same size"	
				return(obj);
			} 
			
			obj@imageList=images;
			obj@dimension=dim(images[[1]]);
			return(obj)
		}
)

# Method to check if moment type is supported
setReplaceMethod (
		f="momentType",
		signature=c("MultiIm"),
		definition = function(obj,value){
			type=tolower(as.character(value));
			if (length(type)==2){		#bivariate polynomials
				if ((sum(type[[1]] == c("cheby","krawt","hahn")) + sum(type[[2]] == c("cheby","krawt","hahn")))==2
					|| (sum(type[[1]] == c("gegen","chebycont","legend")) + sum(type[[2]] == c("gegen","chebycont","legend")))==2){
				obj@momentType=c(type[[1]], type[[2]]);
				} else {
					obj@error[length(obj@error)+1]="This combination of moment types is not supported";
				}
			} else if (length(type)==1){
				if (sum(type==c("gpzm", "fm", "fc", "fr"))==1){
					obj@momentType=type;
				} else if (sum(type==c("krawt", "cheby", "hahn", "chebycont","gegen","legend"))==1){
					obj@momentType=c(type, type);
				} else {
					obj@error[length(obj@error)+1]="The specified moment type is not supported";
				} 
			} else {		
				obj@error[length(obj@error)+1]="The specified moment type is not supported";
			} 
			
			return(obj)
		}
)

# Method to set order
setReplaceMethod(
		f="setOrder",
		signature=c("MultiIm"),
		definition= function(obj, value) {
			if (length(obj@momentType)==2){		#2 moment types
				if (length(value)==2){
					obj@order= as.integer(value);	
				} else if (length(value)==1) {
					obj@order= c(as.integer(value), as.integer(value))
				} else {
					obj@error[length(obj@error)+1] = "Orders specified incorrectly - specify 1 or 2 orders";
				}
			} else if (length(obj@momentType)==1) {
				if (sum(obj@momentType==c("fm", "fc", "fr"))==1){
					if (length(value)==2){
						obj@order= as.integer(value);
					} else {
						obj@error[length(obj@error)+1] = "This moment type requires 2 values: order and repetiton";
					}
				} else if (obj@momentType=="gpzm"){
					if (length(value)==1){
						obj@order= as.integer(value);
					} else {
						obj@error[length(obj@error)+1] = "The moment type requires 1 order";
					}
				} else {
					obj@error[length(obj@error)+1] = "A valid moment type must be specified";
				} 
			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}
			return(obj)
		}		
)

# Method to check necessary Paramsater for each moment type
setReplaceMethod(
		f="setParams",
		signature=c("MultiIm"),
		definition = function(obj,value){
			
			if (length(obj@momentType)==2){			# real orthogonal moments
				if ((length(value)==1) || (length(value)==2)){
					param = list(NULL,NULL);
					for (i in 1:length(value)){
						if (obj@momentType[[i]]=="hahn"){
							if (length(value[[i]])==2){
								param[[i]]= as.numeric(c(value[[i]][[1]], value[[i]][[2]]));
							} else {
								obj@error[length(obj@error)+1] = "Hahn polynomials require 2 parameters: a, c";
							}
						} else if (sum(obj@momentType[[i]]==c("krawt", "gegen"))==1){
							if (length(value[[i]])==1) {
								param[[i]]= as.numeric(value[[i]][[1]])
							} else {
								obj@error[length(obj@error)+1] = paste(obj@momentType[[i]], "requires 1 parameter")
							}
						}
					}
					
					obj@Params= param;
					
				}else {
					obj@error[length(obj@error)+1] = "Parameters incorrectly specified";
					return(obj)
				}
			} else if (length(obj@momentType)==1) {			#complex moments
				if (obj@momentType=="gpzm"){
					a=as.integer(value);
					if (length(a)==1){
						obj@Params=list(a);
					} else {
						obj@error[length(obj@error)+1] = "Generalized Pseudo-Zernike polynomials require 1 parameter";
					}
				} else if (sum(obj@momentType==c("fm", "fc", "fr"))==1){
					obj@Params=list(NULL)
					obj@error[length(obj@error)+1] = "The specified moment type does not require a parameter";
				} else {
					obj@error[length(obj@error)+1] = "A valid moment type must be specified";
				} 
			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}

			return(obj)
		}
)

# Generics for setPoly method
setGeneric (
		name= "setPoly<-",
		def=function(obj,value) standardGeneric("setPoly<-"),
)

# Method to compute polynomials of real orthogonal moment types
setReplaceMethod (
		f="setPoly",
		signature=c("MultiIm"),
		definition = function(obj,value){
			if (length(obj@momentType)==2){
				obj@polynomials= OrthPolynomials(obj@momentType,obj@dimension,obj@order,obj@Params);
			} else if (length(obj@momentType)==1){
				obj@error[length(obj@error)+1] = "This method only available for OrthIm moment types"
			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}
			return(obj)
		}
)


# Method to compute moment
setReplaceMethod (
		f="Moments",
		signature=c("MultiIm"),
		definition = function(obj,value){			
			if (length(obj@momentType)==2){		#real orthogonal moment types
				setPoly(obj)=0;
				xyMat= obj@polynomials;
				xMat=xyMat[[1]];
				yMat=xyMat[[2]];
					
				obj@moments= lapply(obj@imageList,function(z){
								return(t(yMat)%*%z%*%xMat)
							});
			} else if (length(obj@momentType)==1){		#complex moment types
			
					I=array(1, c(obj@dimension[1],obj@dimension[2]))
					# Assume centroid at center of image
					x0 = (obj@dimension[2]+1)/2;
					y0 = (obj@dimension[1]+1)/2;

					radiusTheta=polarXY(I, c(x0,y0));
				
					# Get radius and theta
					radius= radiusTheta[[1]]; 				
					theta= radiusTheta[[2]];

					if (obj@momentType=="gpzm"){
						order= obj@order;
						a = obj@Params[[1]];
						
						temp= GPZMmultiple(obj@imageList,radius,theta,order,a, obj@storePoly);
						obj@moments= temp[[1]]; 
						obj@polynomials= list(temp[[2]]);
						
					} else if (sum(obj@momentType==c("fc","fr", "fm"))) {
						nmax=obj@order[1];
						mmax=obj@order[2];
					
						temp=Cmplxmultiple(obj@momentType,obj@imageList, radius,theta,nmax,mmax, obj@storePoly);
						obj@moments= temp[[1]]; 
						obj@polynomials= list(temp[[2]]);
					} else {
						obj@error[length(obj@error)+1] = "A valid moment type must be specified";
					}

			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}
			
			return(obj)
		}
)

# Method to reconstruct images
setReplaceMethod (
		f="Reconstruct",
		signature= c("MultiIm"),
		definition = function(obj,value){
			#if no maximum order given for reconstruction, use maximum available (not nan or inf)
			order = checkOrder(obj@moments[[1]]);
			#if order defined, check values
			if(!is.null(value)) {
				#check that user-defined order is valid
				if(value[1]>order[1]) {
					obj@error[length(obj@error)+1] = sprintf("maximum order in x direction is %d. Order is being reset to this value.",order[1])
					value[1] = order[1]
				}
				if(obj@momentType[1]!="gpzm") {
					if((length(value)==2) && (value[2]>order[2])) {
						obj@error[length(obj@error)+1] = sprintf("maximum order in y direction is %d. Order is being reset to this value.",order[2])
						value[2] =  order[2]
					} else if ((length(value)!=2)) {
						value[2] =  order[2]
					}
				}
				order=value
			}
			
			if (length(obj@momentType)==2){
				obj@reconstruction= lapply(obj@moments, function(z){
					OrthReconstruct(z, obj@polynomials, obj@dimension, order)
					});
			} else if (length(obj@momentType)==1) {
				# Assume centroid at center of image
				x0 = (obj@dimension[2]+1)/2;
				y0 = (obj@dimension[1]+1)/2;
				I=array(1, c(obj@dimension[1],obj@dimension[2]))
				radiusTheta=polarXY(I, c(x0,y0));
				radius= radiusTheta[[1]]; 			
				theta= radiusTheta[[2]];

				if (obj@momentType=="gpzm"){
					obj@reconstruction= GPZMreconMulti(obj@moments, radius, theta, order[1], obj@Params[[1]]);
				} else if (sum(obj@momentType==c("fc","fr", "fm"))) {
					obj@reconstruction= CmplxreconMulti(obj@momentType, obj@moments, radius, theta, order[1], order[2]);
				} else {
					obj@error[length(obj@error)+1] = "A valid moment type must be specified";
				}
			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}

			return(obj)
		}
)

setReplaceMethod (
		f="Invariant",
		signature= c("MultiIm"),
		definition = function(obj, value){
			#for complex moment types 
			if (length(obj@momentType)==1){
				if (sum(obj@momentType==c("gpzm","fc","fr","fm"))==1){
					if (length(obj@imageList)!= length(obj@moments)){
						Moments(obj)=NULL;
					} 
					obj@invariant= lapply(obj@moments, function(z) {
										z=abs(z);
										return(z/z[1,1])});
				} else {
					obj@error[length(obj@error)+1] = "A valid moment type must be specified";
				}
			} else if (length(obj@momentType)==2) {		#for orthogonal types 
				storeObj= obj;
				transform(obj)= value;
				Moments(obj)= NULL;
				storeObj@invariant= obj@moments;
				obj= storeObj;
			} else {
				obj@error[length(obj@error)+1] = "A valid moment type must be specified";
			}
			
			return(obj)
		}
)

setMethod(
		f="plotPoly",
		signature=c("MultiIm","Numbers"),
		definition = function(obj,order) {
			if (length(dim(order))<2){
				obj@error[length(obj@error)+1] = "Please specify order and repetition to plot"
				return(obj)
			}
			if ((dim(order)[1]<1)||(dim(order)[2]!=2)) {
				obj@error[length(obj@error)+1] = "Please specify order and repetition to plot"
				return(obj)
			}
			
			levels = seq(0,1,.000001);
			g = gray(levels);
			num = dim(order)[1];
			order[,1] = order[,1]+1;
			order[,2] = order[,2]+1;
			
			if (length(obj@momentType)==1){
				if (!(obj@storePoly)) {
					obj@error[length(obj@error)+1] = "obj@storePoly is FALSE. Polynomial may not be stored.";
					return(obj);
				}
 
				if (num<4) {
					par(mfrow=c(1,num),xaxt="n",yaxt="n");
				} else {
					par(mfrow=c(floor(sqrt(num)), ceiling(sqrt(num))),xaxt="n",yaxt="n");
				}
				for (i in 1:num){
					image(asinh(Re(obj@polynomials[[1]][,,order[i,1], order[i,2]])), col=g);
					title(paste("Order", order[i,1]-1, ",Repetition", order[i,2]-1,sep=" "));
				}
			} else {
				if (num<4) {
					par(mfrow=c(1,num),xaxt="n",yaxt="n");
				} else {
					par(mfrow=c(floor(sqrt(num)), ceiling(sqrt(num))),xaxt="n",yaxt="n");
				}				
				for (i in 1:num){
					image(asinh(Re(obj@polynomials[[2]][, order[i,2] +1]%*%t(obj@polynomials[[1]][, order[i,1]+1]))), col=g);
					title(paste("Order in X: ", order[i,1]-1, ",Order in Y: ", order[i,2]-1,sep=" "));
				}
			}
		}
)

setMethod(
		f="plotMoment",
		signature=c("MultiIm"),
		definition = function(obj){
			num = length(obj@moments);
			par(mfrow=c(ceiling(sqrt(num)), num%/%ceiling(sqrt(num))),xaxt="n",yaxt="n");
			for (i in 1:num){
				obj@moments[[i]] -> M;
				if (is.complex(M)) M = abs(M);
				image(asinh(M));
			}
		}
)


setReplaceMethod(
		f="transform",
		signature=c("MultiIm"),
		definition = function(obj, value){
			if (is.null(value) || value<1){
				obj@error[length(obj@error)+1] = "Resolution must be an natural number >0";
				return(obj)	
			}
			
			resolution= as.integer(value);
			
			if (length(value)==1){
				checkSize(obj)= lapply(obj@imageList, function(z){
						return(polarTransform(z, resolution,,calcCentroid(z))[[1]])
						})
			} else {
				scale= as.integer(value[[2]]);
				if (scale<=1){
					obj@error[length(obj@error)+1] = "Scale must be an integer >1";
					return(obj)
				}
				checkSize(obj)= lapply(obj@imageList, function(z){
						return(polarTransform(z, resolution, scale,calcCentroid(z))[[1]])
						})
			}
			
			return(obj)
		}
)












