mww_eval<-function(d,x,filter,LU=NULL){

## mww_eval.m computes the multivariate Wavelet Whittle criterion 
## for the estimation of the long-memory parameter at value d, with exact 
## DWT of Fay et al (2009).
##
## If LU contains different values of lower and upper scales, the minimum
## criterion is returned.
## 
## 	INPUT	d	kx1 long-range dependence parameters
##              x	Data (nxk vector)
##		filter	Wavelet filter
##              LU	Bivariate vector (optional) containing 
##			L, the lowest resolution in wavelet decomposition
##               	U, the maximal resolution in wavelet decomposition
##
##	OUTPUT		Wavelet Whittle criterion
##				
##                                           Achard & Gannaz (2014)
##______________________________________________________________________


if(is.matrix(x)){
	N <- dim(x)[1]
	k <- dim(x)[2]
}else{
	N <- length(x)
	k <- 1
}
x <- as.matrix(x,dim=c(N,k))

## Wavelet decomposition
xwav <- matrix(0,N,k)
for(j in 1:k){
	xx <- x[,j]
	   
	resw <- DWTexact(xx,filter)
	xwav_temp <- resw$dwt
	index <- resw$indmaxband
	Jmax <- resw$Jmax
	xwav[1:index[Jmax],j] <- xwav_temp
}
## we free some memory
new_xwav <- matrix(0,min(index[Jmax],N),k)
if(index[Jmax]<N){
	new_xwav[(1:(index[Jmax])),] <- xwav[(1:(index[Jmax])),]
}
xwav <- new_xwav
index <- c(0,index)

## Wavelet scales
if(is.null(LU)==TRUE){
	LU <- c(1,Jmax)
}
L <- max(LU[1],1)
U <- min(LU[2],Jmax) ## in order to force to be in the correct interval of scales
nscale <- U-L+1

n <- index[U+1]-index[L]

if(k==1){
res<-mww_wav_eval(d,as.vector(xwav),index,LU)
}else{
res<-mww_wav_eval(d,xwav,index,LU)
}

return(res)
}
