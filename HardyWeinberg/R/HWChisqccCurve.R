`HWChisqccCurve` <- function(p,q,chiquant,n,cc=0.5,curvetype="DposUp") {
 switch(curvetype,
 DposUp  = 2*p*q+2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
 DposLow = 2*p*q+2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
 DnegUp  = 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
 DnegLow = 2*p*q-2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n))
}

