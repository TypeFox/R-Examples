soilwater2 <- 
function (x, x0 = 6.653, k0, k1, n) 
	k1 * ( exp(-k0 / x0^n) - exp(-k0 / x^n) )