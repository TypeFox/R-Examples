Hals.snow <-
function(j, Z, Hs, Ht, Hst.ls, b.lag, GP.mx) {
    
	rho <- GP.mx[j, 1]
	reg <- GP.mx[j, 2]
    
    Z.hat <- H.als.b(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, rho=rho, reg=reg, b.lag=b.lag, Hs0=NULL, Ht0=NULL, Hst0.ls=NULL)$Z.hat
   
	return(Z.hat)
}
