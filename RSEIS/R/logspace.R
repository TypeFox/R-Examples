`logspace` <-
function(d1, d2, n=n)
  {
###%LOGSPACE Logarithmically spaced vector.
###%   LOGSPACE(X1, X2) generates a row vector of 50 logarithmically
###%   equally spaced points between decades 10^X1 and 10^X2.  If X2
###%   is pi, then the points are between 10^X1 and pi.
###%
###%   LOGSPACE(X1, X2, N) generates N points.
###%   For N < 2, LOGSPACE returns 10^X2.
###%
###%   See also LINSPACE, :.

###%   Copyright 1984-2002 The MathWorks, Inc. 
###%   $Revision: 5.11 $  $Date: 2002/01/24 06:05:07 $

    if(missing(n)) { n = 50 }

    if(d2 == pi)
      {
        d2 = log10(pi);
      }
    vec = c(d1+(0:(n-2))*(d2-d1)/(floor(n)-1), d2)
    
    
    y = (10)^vec
    
    return(y)
  }

