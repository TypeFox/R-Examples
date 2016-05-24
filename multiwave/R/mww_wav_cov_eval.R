mww_wav_cov_eval<-function(d,xwav,index,psih,grid_K,LU){

## mww_cov_eval_wav.m computes the multivariate wavelet Whittle estimator 
## of the long-run correlation matrix with exact DWT of Fay et al (2009). 
## 
## 	INPUT	 d: kx1 parameter value
##               x: data (nxk vector)
##		 filter: the wavelet filter
##               psih: the Fourier transform of the wavelet mother at
##               values grid_K
##               grid_K: the grid for the approximation of the integral 
##               in K
##               L: the lowest resolution in wavelet decomposition
##               (optional)
##               U: the maximal resolution in wavelet decomposition
##               (optional)
##
##                                           Achard & Gannaz (2014)
##______________________________________________________________________

if(is.matrix(xwav)){
	N <- dim(xwav)[1]
	k <- dim(xwav)[2]
}else{
	N <- length(xwav)
	k <- 1
}
x <- as.matrix(xwav,N,k)

Jmax <- length(index)-1

## Wavelet scales
if(is.null(LU)){
	LU <- c(1,Jmax)
}

  
L <- max(LU[1],1)
U <- min(LU[2],Jmax)
nscale <- U-L+1

n <- index[U+1]-index[L]

sum_xwav <- array(0,dim=c(nscale,k,k))
vect <- matrix(0,nscale,1)
for(f in 1:nscale){
        j <- L+f-1 
        fj <- 1-j
	if(k>1){
        temp <- xwav[(index[j]+1):index[j+1],]%*%diag(2^(fj*d))
			}else{
 temp <- xwav[(index[j]+1):index[j+1]]*2^(fj*d)
				}
        sum_xwav[f,,] <- t(temp)%*%temp
        vect[f] <- fj*(index[j+1]-index[j])
}

    g_temp <- colSums(sum_xwav,1)/n
    g <- matrix(0,k,k)
    g <-  g_temp    
 
## Correction of the phase shift
K <- K_eval(psih,grid_K,d);
if(k==1){
K <- ((exp(-1i*pi*d/2))*K*(exp(1i*pi*d/2)))
}else{
K <- (diag(exp(-1i*pi*d/2))%*%K%*%diag(exp(1i*pi*d/2)))
}
K <- Re(K)

g <- g*(K^-1)

return(g)
}
