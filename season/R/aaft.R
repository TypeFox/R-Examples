# aaft.R
# AAFT (Amplitude Adjusted Fourier Transform), copied from Matlab code
# Jan 2012
# For more details see: D. Kugiumtzis, Surrogate data test for nonlinearity including monotonic transformations, Phys. Rev. E, vol. 62, no. 1, 2000.
# generates surrogate data for a given time series data using the AAFT  
## INPUT:
# - data        : The original time series (column vector)
# - nsur      : The number of surrogate data to be generated 
## OUTPUT:
# - surrogates        : The AAFT surrogate data (matrix of nsur columns)

aaft = function(data, nsur){
n = length(data);
# The following gives the rank order, ixV
ixV=order(data);
rxV=rank(data); # ranks
oxV=data[ixV] # smallest to largest
# ===== AAFT algorithm 
surrogates=matrix(data=0,n,nsur);
for (count in 1:nsur){
   # Rank ordering white noise with respect to data
   rV = rnorm(n);  # random N(0,1)
   orV=rV[order(rV)]; # in order, smallest to largest
   yV = orV[rxV]; # 
   # cor(yV,data,method='spearman') # Not run, should be one
# >>>>> Phase randomisation (Fourier Transform): yV -> yftV 
   if(n%%2==0){n2 = n/2}else{n2 = (n-1)/2}
   tmpV = fft(yV,2*n2); # FFT
   magnV = abs(tmpV); # magnitude
   fiV = Arg(tmpV); # phases
   rfiV = runif(n2-1) * 2 * pi; # random phases
   nfiV = c(0, rfiV, fiV[n2+1], -rev(rfiV));
# New Fourier transformed data with only the phase changed
   tmpV = c(magnV[1:(n2+1)],rev(magnV[2:n2]));
   c.exp=cos(nfiV)+ 1i*sin(nfiV) # complex exponential
   tmpV = tmpV * c.exp; 
# Transform back to time domain;
   yftV=Re(fft(tmpV,inverse=TRUE)); # 3-step AAFT;
# <<<<<
   iyftV = rank(yftV); # Rank ordering data with respect to yftV 
   surrogates[,count] = oxV[iyftV];  # surrogates is the AAFT surrogate of data
} # end of loop
return(surrogates)
} # end of function
