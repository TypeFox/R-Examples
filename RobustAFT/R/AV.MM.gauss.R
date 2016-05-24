"AV.MM.gauss" <- 
function(k0=1.5477, k1=4.6873, sigma=1)
{old <- comval()
 dfcomn(ipsi=4, xk=k1)
 z    <- liepsu(upper=10)
 Q1   <- z$epsi2
 M1   <- z$epsip
 dfcomn2(ipsi=4, xk=k0)
 Beta <- integrate(Chiphi,  -10,10)$value
 Q2   <- integrate(Chi2phi, -10,10)$value-Beta^2
 M2   <- integrate(Chipzphi,-k0,k0)$value
v.lambda <- (sigma^2 * Q1)/M1^2
v.sigma <- (sigma^2 * Q2)/M2^2
dfcomn2(ipsi=old$ipsi, xk=old$xk)
list(V.lambda = v.lambda, V.sigma = v.sigma)}

