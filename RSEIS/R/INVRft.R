`INVRft` <-
function(G,n,tstart,dt)
{
  
  ### function [g,f,t]=INVRft(G,n,tstart,dt);
##%   INVRft           Inverse Fourier Transform, time shifts, sample interval scaling
##% USAGE: [g,f,t]=INVRft(G,n,tstart,dt);
##% Compute the Inverse Fourier Transform of G(f) and order the output to match
##% a time vector that starts at tstart and has a sample interval of dt.
##% n is the length of the time vector and of g.  An N point IFFT is calculated
##% where N is the length of G.  Output is scaled according to the conventions of
##% continuous transforms in Aki and Richards and in J.H. Karl.
##%
##% INPUT: G is a column vector spectrum evaluated at positive and negative
##%        frequencies as defined by MAKEFREQ.
##%        tstart, dt and n define the output time vector as described above.
##%
##% OUTPUT: g is the Inverse Fourier Transform of G.  it is scaled by dt to be
##%        consistent with the continuous transform.  the time shift
##%        theorem has been used to account for time not starting at t=0.
##%        f and t are the time and frequency vectors for g and G.
##%        the lengths of g and t are n.
##%
##% See also FRWDFT, and MAKEFREQ.

##% adapted and translated from K. Creager  kcc@geophys.washington.edu   12/30/93
########  modified by J. M. Lees 10/20/2007


N=length(G);                    ##% length of input time-domain vector
## t=tstart(1)+[0:dt:(n-1)*dt]';   ##% time vector
t=tstart[1]+seq(from=0,by=dt,to=(n-1)*dt)
f=makefreq(N,dt)
i = complex(real=0, imaginary=1)
g=(1/(N*dt))*fft(G*exp(2*pi*i*tstart*f), inverse = TRUE); ##% inverse fourier transform with time shift
g=g[1:N];                             ##%truncate time vector to n points.

return(list(g=g, f=f, t=t))


}

