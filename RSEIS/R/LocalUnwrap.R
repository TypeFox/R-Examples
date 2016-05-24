`LocalUnwrap` <-
function(p,cutoff=cutoff)
{
###%LocalUnwrap   Unwraps column vector of phase values.
if(missing(cutoff)) { cutoff = pi;  }
m = length(p); 

###% Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
###% by constraining it to the range [-pi,pi]
dp = diff(p);
##% Incremental phase variations
dps = (  RPMG::fmod(dp+pi, 2*pi) ) - pi;      ##% Equivalent phase variations in [-pi,pi)




dps[dps==-pi & dp>0] = pi;     ##% Preserve variation sign for pi vs. -pi


dpcorr = dps - dp;              ##% Incremental phase corrections

dpcorr[abs(dp)<cutoff] = 0;   ##% Ignore correction when incr. variation is < CUTOFF

##%% Integrate corrections and add to P to produce smoothed phase values
p[2:m] = p[2:m] + cumsum(dpcorr);

return(p)

}

