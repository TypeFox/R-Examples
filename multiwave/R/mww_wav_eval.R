mww_wav_eval<-function(d,xwav,index,LU=NULL){
## multi_wav_exact.m computes the multivariate Wavelet Whittle criterion 
## for the estimation of the long-memory parameter at value d, with exact 
## DWT of Fay et al (2009).
##
## If LU contains different values of lower and upper scales, the minimum
## criterion is returned.
## 
## 		INPUT	d: kx1 parameter value
##               	x: data (nxk vector)
##			xwav: the wavelet coefficients
##               	LU: the first and maximal resolution level used for 
##               estimation in wavelet decomposition (optional)
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
Jmax <- length(index)-1

## Wavelet scales
	if(is.null(LU)){
		LU <- c(1,Jmax)
	}
			
    L <- max(LU[1],1)
    U <- min(LU[2],Jmax) ## in order to force to be in the correct intreval of scales
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
        	sum_xwav[f,,] <- t(temp)%*%temp;
        	vect[f] <- fj*(index[j+1]-index[j]);
	}

    g_temp <- colSums(sum_xwav,1)/n
    g <- matrix(0,k,k)
    g <-  g_temp
    g  <- Re(g);    
 
    ## criterion
    r <- log(det(g)) - 2*sum(d)*log(2)*sum(vect)/n

return(r)
}
