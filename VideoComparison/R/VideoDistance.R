########################################################################
# ph_hammingdistance2 is translated into R from phash.cpp (http://www.phash.org/) to avoid installing all the library and its dependences.
# ph_crosscorr uses a different approach from pHash to obtain the same results
# We would like to thank Evan Klinger & David Starkweather for their library
# These functions are not used because not all the phash distance methods were translated into R and the c code is still required 
ph_hammingdistance2<-function(hashA,hashB){ #hashA & hashB are integer vectors and their elements should be in 0:255 range
	##Jose M Perez Ramos: Translated from phash.cpp and using R's vectorized functions

	if(length(hashA)!=length(hashB))
		return(-1.0)

	if(length(hashA)<=0)
		return(-1.0)
	
	ph_bitcount8<-function(val){
		num = integer(length=length(val));

		repeat{
			greaterthan0 = val>0
			if(!any(greaterthan0))
				break

			num[greaterthan0] = num[greaterthan0] +1L
		  	val = bitwAnd(val,val-1)
		}
		return(sum(num))
	}

	xored = bitwAnd(bitwXor(hashA,hashB),0xFF) #force hashA & hashB are in 0:255 range
	return(ph_bitcount8(xored)/(8*length(hashA)))
}
ph_crosscorr<-function(x_coeffs,y_coeffs){ #x_coeffs and y_coeffs are integer vectors and their elements should be in 0:255 range
    ##Jose M Perez Ramos: Using FFT to implement the xcorr function. FFT requires no padding because it is a 'cyclic' xcorr

    if(length(x_coeffs)!=length(y_coeffs))
    	return(-1)

    x = x_coeffs - mean(x_coeffs)
    y = y_coeffs - mean(y_coeffs)
    res = Re((fft(fft(x)*Conj(fft(y)), inverse = TRUE)/length(y))/sqrt(sum(y^2)*sum(x^2))) # cyclic xcorr
    return(max(res))
}
########################################################################

VideoDistance<-function(hh,h2) {
	if ( ! is.list(hh) || ! is.list(h2) ) {
      rs<-list(dct=-1,mw=-1,str=-1,rd=-1,errlevel=1,errtxt="No hash provided")
		return(rs)
   }
	res<-.Call("VideoDistance",hh,h2)
   rs<-list(dct=as.integer(res[1]), str=as.integer(res[4]),mw=res[2],
			rd=res[3], errlevel=0,errtxt="")
	return(rs);
}
