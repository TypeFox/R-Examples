#EZRL.R
#source for functions: EZR, EZL, Erf, Erfc

EZR <- function(zmu, zsig, c){
  (zmu - (sqrt(2/pi)*zsig)/exp((c - zmu)**2/(2.*zsig**2)) + 
     zmu*Erf((c - zmu)/(sqrt(2)*zsig)))/Erfc((-c + zmu)/(sqrt(2)*zsig))
}

EZL <- function(zmu, zsig, c){
(sqrt(2/pi)*zsig + 
   exp((c - zmu)**2/(2.*zsig**2))*zmu*Erfc((c - zmu)/(sqrt(2)*zsig)))/
  (exp((c - zmu)**2/(2.*zsig**2))*(1 + Erf((-c + zmu)/(sqrt(2)*zsig))))
}

Erf <- function(z) -1 + 2*pnorm(sqrt(2)*z)
Erfc <- function(z) 2*pnorm(-(sqrt(2)*z))