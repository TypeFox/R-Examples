`radSV` <-
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

A1 = sin(lam)*cos(2*del)*cos(2*ichi)*sin(phidif)
A2 = cos(lam)*cos(del)*cos(2*ichi)*cos(phidif)
A3 = 0.5*cos(lam)*sin(del)*sin(2*ichi)*sin(2*phidif)
A4 = 0.5*sin(lam)*sin(2*del)*sin(2*ichi)*(1+(sin(phidif)*sin(phidif)))


FSV = A1 -A2 +A3 -A4


#  FSV = (sin(lam)*cos(2*del)*cos(2*ichi)*sin(phidif)-cos(lam)*cos(del)*cos(2*ichi)*cos(phidif)+ 0.5*cos(lam)*sin(del)*sin(2*ichi)*sin(2*phidif)-0.5*sin(lam)*sin(2*del)*sin(2*ichi)*(1+(sin(phidif)*sin(phidif))))

return(FSV)

}

