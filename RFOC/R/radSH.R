`radSH` <-
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
FSH = cos(lam)*cos(del)*cos(ichi)*sin(phidif)+cos(lam)*sin(del)*sin(ichi)*cos(2*phidif)+ sin(lam)*cos(2*del)*cos(ichi)*cos(phidif)-0.5*sin(lam)*sin(2*del)*sin(ichi)*sin(2*phidif)

return(FSH)

}

