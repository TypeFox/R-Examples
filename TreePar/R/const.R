const <-
function(i,t,l,mu,rho) {
	if (i==1) {
		ci = 1-rho[1]
		} else {
		ci = 1-rho[i] + rho[i] * q2(i-1,t[i],t,l,mu,rho)	}
	ci
	}

