`radP` <-
function( del, phiS, lam, ichi, phi)
{

#  convert all angles to radians
DEG2RAD = pi/180

lam = DEG2RAD*lam
del = DEG2RAD*del
phiS = DEG2RAD*phiS
ichi = DEG2RAD*ichi
phi  = DEG2RAD*phi

phidif  = phi - phiS
Fp = cos(lam)*sin(del)*sin(ichi)*sin(ichi)*sin(2*phidif)-cos(lam)*cos(del)*sin(2*ichi)*cos(phidif)+ sin(lam)*sin(2*del)*(cos(ichi)*cos(ichi)-sin(ichi)*sin(ichi)*sin(phidif)*sin(phidif))+sin(lam)*cos(2*del)*sin(2*ichi)*sin(phidif)

return(Fp)

}

