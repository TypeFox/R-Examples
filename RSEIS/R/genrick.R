genrick <-
function(freq,dt,nw)
{
#
# RICKER WAVELET GENERATOR
# Written by Leonard Lisapaly (leonardl@fisika.ui.ac.id)
  ##  converted to R  j. lees Wed Sep  1 09:39:21 2004
#
# INPUTS
# freq = wavelet dominant frequency [Hz]
# dt   = sampling interval [sec]
# nw   = length of wavelet [odd number]
#
# USAGE
# If you want to obtain a 35 samples Ricker wavelet 
# with dominant frequency of 25 Hz and sampling interval
# of 0.002 sec, you should type :
#      w = genrick(25,0.002,35)

a  = freq*sqrt(pi)/2;
nc = (nw+1)/2;
tc = (nc-1)*dt;
t  = seq(from=0, length=nw-1 )*dt;
b  = pi*freq*(t-tc);
w  = a*(1-2*b^2)*exp(-b^2);
return(w)

}

