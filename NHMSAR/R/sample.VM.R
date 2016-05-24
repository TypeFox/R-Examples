sample.VM <-
function(a,b) {
	
	tau = 1+sqrt(1+4*b^2)
	rho = (tau-sqrt(2*tau))/(2*b)
	r = (1+rho*rho)/(2*rho)
	while (1) {
		u1 = runif(1) 
		z = cos(pi*u1)
		f = (1+r*z)/(r+z)
		c = b*(r-f)
		u2 = runif(1)
		if (u2<c*(2-c)) {break}
		if (c<=log(c/u2)+1) {break}
	}
	u3 = runif(1)
	if (u3<.5) {x = a-acos(f)}
	else {x = a+acos(f)}
	return(x)
	
}
