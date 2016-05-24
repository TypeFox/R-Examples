# General use functions, these are used by image class
# 
# Author: Allison Irvine, Tan Dang
###############################################################################


#continuous complex moments 
#gpzm		Generalized Pseudo-Zernike
#fm			fourier mellin
#fc			fourier chebyshev
#fr			radial fourier

#orthogonal moments - discrete or continuous
####### can choose a combination of 2 from the same group: 
#DISCRETE:
#"cheby"		Chebyshev
#"krawt"		Krawtchouk
#"hahn"			Hahn

#CONTINUOUS:
#"gegen"		Gegenbauer
#"chebycont"	ChebyshevCont	(continuous chebyshev)
#"legend"		Legendre

setClassUnion("Numbers", c("numeric", "array", "matrix"))

#standalone method calculates moments and returns object
#requires arguments: image,momentType,order,params
setGeneric (
		name = "momentObj",
		def=function(I,type,order,paramsX,paramsY) {standardGeneric("momentObj")},
)

setMethod(
		f="momentObj",
		signature=c("Numbers","character","Numbers"),
		definition = function(I,type,order,paramsX,paramsY) {
			
			if (missing(paramsX)) paramsX= NULL;
			if (missing(paramsY)) paramsY= NULL;
			
			#allow for pairs of orthogonal moments
			if(length(type)>1){
				if(((sum(type[1]==c("cheby","krawt","hahn"))+sum(type[2]==c("cheby","krawt","hahn")))==2) 
					||	((sum(type[1]==c("gegen","chebycont","legend"))+sum(type[2]==c("gegen","chebycont","legend")))==2))	{
					#compute bivariate orthogonal moments
					obj = new("OrthIm",img=I);
					momentType(obj) = type;
					setOrder(obj) = order;
					setParams(obj) = list(paramsX, paramsY);
					Moments(obj) = NULL;
				}
			}else if((sum(type==c("cheby","krawt","hahn","gegen","chebycont","legend")))==1){
					#compute orthogonal moments
					obj = new("OrthIm",img=I);
					momentType(obj) = type;
					setOrder(obj) = order;
					setParams(obj) = list(paramsX,paramsX);
					Moments(obj) = NULL;
			} else if(sum(type==c("gpzm","fm","fc","fr"))==1){	
					#complex moments
					obj = new("CmplxIm",img=I)
					momentType(obj) = type
					setOrder(obj) = order;
					setParams(obj) = paramsX;
					Moments(obj) = NULL;
			}else{
				return("A valid moment type must be specified")
			}
			
			return(obj)
		}
)

#defined for a list of images
setMethod(
		f="momentObj",
		signature=c("list","character","Numbers"),
		definition = function(I,type,order,paramsX,paramsY) {
			
			if (missing(paramsX)) paramsX= NULL;
			if (missing(paramsY)) paramsY= NULL;
			
			#allow for pairs of orthogonal moments
			if(length(type)==2){
				if(((sum(type[1]==c("cheby","krawt","hahn"))+sum(type[2]==c("cheby","krawt","hahn")))==2) 
						||	((sum(type[1]==c("gegen","chebycont","legend"))+sum(type[2]==c("gegen","chebycont","legend")))==2))	{
					#compute bivariate orthogonal moments
					obj = new("MultiIm",I);
					momentType(obj) = type;
					setOrder(obj) = order
					setParams(obj) = list(paramsX, paramsY);
					Moments(obj) = NULL;
				}
			}else if(((sum(type==c("cheby","krawt","hahn")))==1) || (sum(type==c("gegen","chebycont","legend"))==1)){
				#compute orthogonal moments
				obj = new("MultiIm",I);
				momentType(obj) = type;
				setOrder(obj) = order;
				if (sum(type==c("krawt","hahn","gegen"))==1){
					setParams(obj) = list(paramsX,paramsX);
				}
				Moments(obj) = NULL;
			} else if(sum(type==c("gpzm","fm","fc","fr"))==1){	
				#complex moments
				obj = new("MultiIm",I)
				momentType(obj) = type;
				setOrder(obj) = order;
				setParams(obj) = paramsX;
				Moments(obj) = NULL;
			}else{
				return("Error")
			}
			
			return(obj)
		}
)


#display grayscale image
setGeneric (
		name= "displayImg",
		def=function(img) standardGeneric("displayImg"),
)

setMethod(
		f="displayImg",
		signature=c("Numbers"),
		definition = function(img) {
			#if image is not grayscale, convert to grayscale
			if(length(dim(img))>2) {
				img = rowSums(img, dims=2)/3
			}
			if(length(dim(img))==2) {
				levels = seq(0,1,.0000001);
				g = gray(levels);
				#rotate image so that it appears aligned
				img = rotate270(img);
				#perform histogram equalization on displayed image
				img <- histeq(img);
				par(mfrow = c(1,1))
				image(img,col=g,axes=FALSE);
			} else {
				return("problem with image format")
			}
		}
)

#to display a list of grayscale images
setMethod(
		f="displayImg",
		signature=c("list"),
		definition = function(img){
			err=FALSE;
			levels = seq(0,1,.0000001);
			g = gray(levels);
			nImg = length(img);
			if (nImg < 4){
				par(mfrow = c(1,nImg),xaxt="n",yaxt="n")
			}else {
				par(mfrow=c(floor(sqrt(nImg)),ceiling(sqrt(nImg))),xaxt="n",yaxt="n");
			}
			for (i in 1:nImg){
				#if image is not grayscale, convert to grayscale
				if(length(dim(img[[i]]))>2) {
					img[[i]] = rowSums(img[[i]], dims=2)/3
				}
				if(length(dim(img[[i]]))==2) {
					temp = img[[i]];
					#perform histogram equalization on displayed image
					temp = histeq(temp)
					#rotate image so that it appears aligned
					image(rotate270(temp),col=g);
				} else {
					err=TRUE;
				}
			}
			
			if (err) {
				return("problem with image format")
			}
		}
)

#rotate image so that it appears aligned
setGeneric (
		name = "rotate270",
		def=function(img) {standardGeneric("rotate270")},
)

setMethod(
		f="rotate270",
		signature=c("Numbers"),
		definition = function(img) {
			im <- as.data.frame(t(img));
			im <- rev(im);
			im <- as.matrix(im);
			return(im)
		}
)
	
	
#performs histogram equalization on image matrix
setGeneric (
		name = "histeq",
		def=function(I) {standardGeneric("histeq")},
)

setMethod(
		f="histeq",
		signature=c("Numbers"),
		definition = function(I) {
			I = (I-min(I))/(max(I)-min(I))
			I = round(I*255);
			
			G =256;	  
			H =array(0:255,256);		
			T =array(0,256);
			
			H = apply(H,1, function(z){ sum(I==z) });
			
			for (i in 2:length(H)){
				H[i]= H[i-1]+H[i]	
			}
			
			T = H*(G-1)/length(I);
			
			for (i in 1:length(I)){
				I[i]=T[I[i]+1]
			}
			
			return(I)
		}
)

# zero out unused moment
zeroMoment<- function(M){
	M=flip(M); 
	M[upper.tri(M,FALSE)]=0; 
	flip(M)
}

#flip the matrix horizontally
flip<- function(M){
	for (i in 1:floor(dim(M)[1]/2)){
		temp=M[i,]
		M[i,]=M[dim(M)[1]-i+1,]
		M[dim(M)[1]-i+1,]=temp
	}
	M
}


#map 1 dimension of image into [-1,1]x[-1,1]
mapCenter <- function(N,centerX,Invariant){
	if (!Invariant){
		x= 0:(N-1);
		return(-1+(x+1/2)*(2/N));
	}
	#x= 0:(N-1);
	#x= (-1+(x+1/2)*(2/N))/(sqrt(Invariant))
	
	#return(x)
	x=0:(N-1);
	
	#max radius on x direction
	maxX= max(abs(1-centerX),abs(N-centerX))
	
	# map x into [-1,1]
	x = ((x-centerX)/maxX)/sqrt(Invariant)
	
	return(x)
}


#function definition for centroid calculation
setGeneric (
		name= "calcCentroid",
		def=function(I) {standardGeneric("calcCentroid")},
)

#this is a general-purpose version of function
setMethod(
		f="calcCentroid",
		signature=c("matrix"),
		definition = function(I) {
			m<- mat.or.vec(3,1);
			NY <- dim(I)[1];
			NX <- dim(I)[2];
			
			#m(0,0)
			m[1] = sum(I);
			#m(0,1)
			m[2] = sum(t(I)%*%(1:NY));
			#m(1,0)
			m[3] = sum(I%*%(1:NX));
			
			#calculate x0 and y0
			y0 = m[2]/m[1];
			x0 = m[3]/m[1];
			
			return(c(x0,y0))
		}
)



#display polynomials of orthogonal moments
setGeneric (
		name = "demoPoly",
		def=function(order,type,N,params) {standardGeneric("demoPoly")},
)

setMethod(
		f="demoPoly",
		signature=c("Numbers","character","numeric"),
		definition = function(order,type,N,params) {
			if (missing(params)) params= NULL;
			
			if (sum(type==c("gpzm", "fm", "fr", "fc"))==1) {	#if complex
				img= list(mat.or.vec(N, N));
				if(length(dim(order))<2) {
					return(NULL)
				}
				#gpzm only requires order, repetition is constrained by order
				if (type=="gpzm") {
					maxOrder= max(order)
				}else {				
					maxOrder= c(max(order[,1]), max(order[,2]));
				}
				
				n = dim(order)[[1]];
				
				obj= new("MultiIm", img)
				momentType(obj)= type;
				obj@storePoly= TRUE;
				setOrder(obj) = maxOrder;
				setParams(obj)= params;
				Moments(obj) = NULL;

				levels = seq(0,1,.000001);
				g = gray(levels);
				if (n<4) {
					par(mfrow=c(1,n),xaxt="n",yaxt="n");
				} else {
					par(mfrow=c(floor(sqrt(n)), ceiling(sqrt(n))),xaxt="n",yaxt="n");
				}
				for (i in 1:n){
					image(asinh(Re(obj@polynomials[[1]][,,(order[i,1]+1), (order[i,2]+1)])), col=g);
					title(paste("Order", order[i,1], "and Repetition", order[i,2],sep=" "));
				}
				
			} else { 	#if orthogonal, calculate polynomials directly
				P=switch(type, 
					krawt=Krawtchouk(max(order),N,params),
					cheby=Chebyshev(max(order),N),
					hahn=DualHahn(max(order),N,params),
					gegen = Gegenbauer(max(order),N,params),
					chebycont = ChebyshevCont(max(order),N),
					legend = Legendre(max(order),N),
				)
			
			
				typeName = ""
				#plot up to 20 colors
				colorVec = c(31,259,91,142,10,12,17,24,26,32,33,43,47,51,54,60,62,73,75,79) 
				
				#get moment type full name
				typeName = switch(type,
						gpzm = "Generalized Pseudo-Zernike",
						fc = "Fourier Chebyshev",
						fr = "Radial Fourier",
						fm = "Fourier Mellin",
						cheby = "Discrete Chebyshev",
						chebycont = "Continuous Chebyshev",
						legend = "Legendre",
						gegen = "Gegenbauer",
						krawt = "Krawtchouk",
						hahn = "Dual Hahn"
				)
					
				#plot all orders
				if (!is.null(order)) {
					X = P[,order[1]];
					plot(1:length(X),X,type="n",xlab="X coordinate",ylab="polynomial(X)")
					#title(main = sprintf("%s Polynomials, order %d to %d",typeName,order[1],order[length(order)]))
					title(main = paste("order",order[1],"to",order[length(order)], sep=" "))
					h = 1;
					for(j in order) {
						X = P[,j];
						lines(1:length(X),X,col=colors()[colorVec[h%%20]])
						h=h+1
					}
				}		
			}#end orthogonal case
		}
)

#check maximum available order of calculated moments
setGeneric (
		name= "checkOrder",
		def=function(moments) standardGeneric("checkOrder"),
)

setMethod(
		f="checkOrder",
		signature=c("matrix"),
		definition = function(moments) {
			NX = dim(moments)[2];
			NY = dim(moments)[1];
			
			valid = !((is.infinite(moments)) | (is.nan(moments)));
			#find largest rectangular block that contains valid values
			
			orderX = max(which(valid[1,]));
			orderY = max(which(valid[,1]));
			
			if(sum(!valid[1:orderY,1:orderX]) > 0) {
				#maximum available order (not containing inf or NaN)
				N = min(NX,NY);
				for (i in 2:N) {
					isInf = apply(moments[1:i,1:i], 1, function(x) is.infinite(x));
					isNan = apply(moments[1:i,1:i], 1, function(x) is.nan(x));
					if ( (sum(sum(isInf)) > 0) | (sum(sum(isNan)) > 0) ) {
						maxOrder = i-1;
						orderX = orderY = maxOrder;
						break;
					} else if(i==N) {
						maxOrder=i;
					}
					
				}
			}
			
			return(c(orderX-1,orderY-1))
		}
)


#pixel coordinates in polar representation 
setGeneric (
		name= "polarXY",
		def=function(I,centroid) {standardGeneric("polarXY")},
)

setMethod(
		f="polarXY",
		signature=c("matrix","numeric"),
		definition = function(I,centroid) {
			x0 = centroid[1];
			y0 = centroid[2];
			
			#get maximum radius
			maxRadius = calcMaxRadius(I, centroid)
			
			#calculate relative x and y coordinates
			y = (((1:dim(I)[1])-y0)/maxRadius) * -1;
			x = ((1:dim(I)[2])-x0)/maxRadius;
			#calculate radius r and angle theta
			theta = atan(y%*%t(1/x));
			r = sqrt(apply(as.matrix(x^2), 1, function(z) z+(y^2)));
			
			#adding pi to quadrants 2 and 3 of theta
			temp = which(x < 0);
			theta[,temp] = theta[,temp] + pi;
			
			result = list(r,theta);
			names(result)=c("radius","theta")
			return(result)
		}
)



#function definition for maximum radius calculation
setGeneric (
		name= "calcMaxRadius",
		def=function(I, center) {standardGeneric("calcMaxRadius")},
)

setMethod(
		f="calcMaxRadius",
		signature=c("matrix"),
		definition = function(I, center) {
			if (missing(center)){
				d = dim(I);
			
				maxRadius= sqrt(sum(I))/2*sqrt(d[1]/d[2]+d[2]/d[1]);
			
				return(maxRadius);
			}
			x0 = center[1];
			y0 = center[2];
			NY = dim(I)[1];
			NX = dim(I)[2];
			
			radius=c(1:4);
			# max radius is the distance of one of the 4 corners of the image to the centroid
			radius[1]= sqrt((1 - x0)^2 + (1 - y0)^2);
			radius[2]= sqrt((1 - x0)^2 + (NY - y0)^2);
			radius[3]= sqrt((NX - x0)^2 + (1- y0)^2);
			radius[4]= sqrt((NX - x0)^2 + (NY - y0)^2);
			
			return(max(radius))

		}
)
#create plot of image intensities as radius vs angle using center of image as default
setGeneric (
		name= "polarTransform",
		def=function(I,resolution,scale,center) standardGeneric("polarTransform"),
)

setMethod(
		f="polarTransform",
		signature=c("matrix","numeric"),
		definition = function(I,resolution,scale,center) {
		
			if (missing(center)) {
				#automatically set center to center of image, maybe we can make this an option	
				center = c(round(dim(I)[2]/2),round(dim(I)[1]/2))
				center= round(center);
			} 
			
			#radius of largest circle with origin at centroid which will fit within image boundaries
			maxR = 	min(abs(dim(I)[1]-center[2]),center[1],abs(dim(I)[2]-center[1]),center[2]);
			maxR = round(maxR)-1;
	
			maxRold=maxR;
			if (!missing(scale)) {
				maxR=scale;
			} else {
				scale=NULL;
			}
		
			#find pixels at the given polar coordinates
			X = mat.or.vec(maxR,maxR*resolution);
			Y = mat.or.vec(maxR,maxR*resolution);
		
			#find principal axis, and shift
			pAxis= getPAxis(I);
		
			#get coordinates relative to centroid
			for (r in 1:maxR) {
				theta= seq(0,by=2*pi/(resolution*r),length.out=(resolution*r))
				X[r,1:(r*resolution)] = r*cos(theta-pAxis)*(maxRold/maxR);
				Y[r,1:(r*resolution)] = r*sin(theta-pAxis)*(maxRold/maxR);
			}
		
			#get actual image coordinates
			X2 = round(X) + center[1] + 1;
			Y2 = round(Y) + center[2] + 1;
		
			xLim = dim(X2)[2];
			yLim = dim(X2)[1];
			PI = mat.or.vec(xLim,yLim);
		
			for (y in 1:yLim){
				x= as.array(1:(y*resolution))
				PI[x,y]= apply(x, 1, function(z){
							return(I[Y2[y,z], X2[y,z]])
						})
			}
			
			PI= t(PI)
		
			#return information neccessary for reverse transform
			result = list(PI,pAxis,resolution,center,scale);
			names(result) = c("PI","pAxis","resolution","center","scale")
			return(result)
	}
)


#create plot of image intensities as radius vs angle using center of image as default
setGeneric (
		name= "revPolar",
		def=function(d,params) standardGeneric("revPolar"),
)

setMethod(
		f="revPolar",
		signature=c("numeric","list"),
		definition = function(d,params) {
			PI = params[[1]]
			pAxis = params[[2]]
			resolution = params[[3]]
			center = params[[4]]
			scale = params[[5]]
			
			maxR = 	min(abs(d[2]-center[1]),center[1],abs(d[1]-center[2]),center[2]);
			
			maxR= round(maxR)-1;
			
			maxRold=maxR;
			if (!missing(scale) || !is.null(scale)) {
				maxR=scale;
			} else {
				scale=NULL;
			}
			
			#find pixels at the given polar coordinates
			X = mat.or.vec(maxR,maxR*resolution);
			Y = mat.or.vec(maxR,maxR*resolution);
			
			for (r in 1:maxR) {
				theta= seq(0,by=2*pi/(resolution*r),length.out=(resolution*r))
				X[r,1:(r*resolution)] = r*cos(theta-pAxis)*(maxRold/maxR);
				Y[r,1:(r*resolution)] = r*sin(theta-pAxis)*(maxRold/maxR);
			}
			
			X2 = round(X) + center[1] + 1;
			Y2 = round(Y) + center[2] + 1;
			
			xLim = dim(PI)[2];
			yLim = dim(PI)[1];
			I2 = mat.or.vec(max(Y2),max(X2))
			
			for (x in 1:xLim) {
				for(y in 1:yLim) {
					I2[Y2[y,x],X2[y,x]] = PI[y,x]
				}
			}
			
			I2
				
		}
)

# These 3 functions are used to normalize image at polarTransfrom
###################
# Compute (central) Geometric Moment
#centralized - option, calculate central moments
geoMoment<- function(n,m,I,centralized){
	
	if (centralized){
		center=calcCentroid(I);
		x= 1:dim(I)[2]-center[1];
		y= 1:dim(I)[1]-center[2];
	} else {
		x= 1:dim(I)[2];
		y= 1:dim(I)[1];
	}
	return(as.numeric(t(y^m)%*%I%*%x^n))
}

# Compute central geometric moment that is invariant to rotation, scaling, and translation
#n,m is the order
#theta is the principal axis
vlm<- function(n,m,theta,I){
	center= calcCentroid(I);
	r= (n+m)/2+1;
	
	#x=1:dim(I)[2];
	y=1:dim(I)[1];
	
	v=0;
	for (x in 1:dim(I)[2]){
		v= v+ sum(((x-center[1])*cos(theta)+(y-center[2])*sin(theta))^n*
					((y-center[2])*cos(theta)-(x-center[1])*sin(theta))^m*I[,x]);
	}
	v=v*sum(I)^(-r);
	
	v
}

#Determine the principal axis
getPAxis<-function(I) {
	#calculate first principal axis, if conditions are not satisfied then add pi/2 until conditions are satisfied 
	theta= 1/2*atan(2*geoMoment(1,1,I, TRUE)/(geoMoment(2,0,I, TRUE)-geoMoment(0,2,I, TRUE)));
	n=1;
	while(!(vlm(2,0,theta,I)>vlm(0,2,theta,I) && vlm(3,0,theta,I)>0)){
		theta= theta+pi/2;
		n=n+1;
		if (n==6) {
			print("Unable to find principal axis")
			return(theta);
		}
	}
	return(theta)
}


