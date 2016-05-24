dsim <-
function (x, mu, sig) {
   	dxmu <- (x-mu)^2/(x*(1-x)*mu^2*(1-mu)^2)
   	return(1/sqrt(2*pi*sig^2*(x*(1-x))^3)*exp(-1/2/sig^2*dxmu))
}
