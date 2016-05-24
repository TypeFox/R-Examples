scaling_filter <- function(family,parameter){

## Computes the filter coefficients of the Daubechies wavelet family
##
## 	INPUT	family  	Wavelet family, 'Haar' or 'Daubechies'
##        	parameter   	Order of the Daubechies wavelet 
##				(equal to twice the number of vanishing moments)
##
## 	OUTPUT  h	Vector of scaling filter coefficients
##         	M 	Number of vanishing moments
##		alpha	Fourier decay exponent
##
##					based on the paper of Fay, Moulines, Roueff, Taqqu, 2009
##                                      Achard & Gannaz (2014)
##________________________________________________________________________________________________

    
choosen <- FALSE
if(family=='Haar'){
    h <- c(1,1)/sqrt(2)
    M <- 1
    alpha <- 1
    choosen <- TRUE
}
else{
    if(family=='Daubechies'){
	if(parameter==2){
		h <- c(1,1)/sqrt(2)
        	M <- 1
        	alpha <- 1
  		choosen <- TRUE
      	}
    	if(parameter==4){
		h <- c(4.829629131445342e-01,
	               8.365163037378079e-01,
	               2.241438680420134e-01,
	       	      -1.294095225512604e-01);
        	M <- 2
        	alpha <- 1.3390
  		choosen <- TRUE
      	}
     	if(parameter==6){
	  	h <- c(3.326705529500830e-01,
	               8.068915093110932e-01,
	       	       4.598775021184915e-01,
	       	      -1.350110200102552e-01,
	       	      -8.544127388202682e-02,
	               3.522629188570960e-02)        
        	M <- 3
        	alpha <- 1.6360
  		choosen <- TRUE
      	}
      	if(parameter==8){
		h <-c(2.303778133088966e-01,
	      	      7.148465705529161e-01,
	      	      6.308807679298595e-01,
	      	     -2.798376941685971e-02,
	      	     -1.870348117190937e-01,
	      	      3.084138183556036e-02,
	      	      3.288301166688513e-02,
	      	     -1.059740178506902e-02) 
	  	M <- 4
		alpha <- 1.9125
  		choosen <- TRUE
      	}
      	if(parameter==10){
		h <- c(1.601023979741941e-01,
	       	       6.038292697971941e-01,
	               7.243085284377780e-01,
	               1.384281459013204e-01,
	              -2.422948870663870e-01,
	              -3.224486958464181e-02,
	               7.757149384004471e-02,
	              -6.241490212798629e-03,
	              -1.258075199908210e-02,
	               3.335725285473789e-03)
        	M <- 5
        	alpha <- 2.1766
		choosen <- TRUE
      	}
      	if(parameter==12){
		h <- c(1.115407433501186e-01,
	       	       4.946238903984910e-01,
	               7.511339080211437e-01,
	       	       3.152503517091979e-01,
	       	      -2.262646939654832e-01,
	              -1.297668675672881e-01,
	               9.750160558731831e-02,
	       	       2.752286553029524e-02,
	              -3.158203931749558e-02,
	               5.538422011600732e-04,
	               4.777257510945812e-03,
	              -1.077301085308614e-03)
        	M <- 6
        	alpha <- 2.4322
  		choosen <- TRUE
      	}
      	if(parameter==14){
		h <-c(7.785205408501074e-02,
	      	      3.965393194819243e-01,
	              7.291320908462448e-01,
	              4.697822874051911e-01,
	             -1.439060039285805e-01,
	             -2.240361849938853e-01,
	              7.130921926683061e-02,
	              8.061260915108569e-02,
	             -3.802993693501259e-02,
	             -1.657454163066412e-02,
	              1.255099855610169e-02,
	              4.295779729216693e-04,
	             -1.801640704047523e-03,
	              3.537137999745342e-04)
        	M <- 7
        	alpha <- 2.6817
  		choosen <- TRUE
      }
      if(parameter==16){
		h <- c(5.441584224313754e-02,
	       	       3.128715909144829e-01,
	       	       6.756307362976408e-01,
	       	       5.853546836543904e-01,
	       	      -1.582910525661015e-02,
	       	      -2.840155429619162e-01,
	       	       4.724845738282760e-04,
	      	       1.287474266205178e-01,
	       	      -1.736930100185194e-02,
	       	      -4.408825393083552e-02,
	       	       1.398102791740789e-02,
	       	       8.746094047410769e-03,
	       	      -4.870352993456092e-03,
	       	      -3.917403733773278e-04,
	       	       6.754494064510721e-04,
	       	      -1.174767841248659e-04)
        	M <- 8
        	alpha <- 2.9265
  		choosen <- TRUE
     }
     

     if(!choosen){ stop('Unavailable order (Daubechies wavelet)') }

     h <- h/norm(as.matrix(h),'F')
  }
}

if(!choosen){ stop('Unknown family of scaling functions') }

list(h=h,M=M,alpha=alpha)
  
}

