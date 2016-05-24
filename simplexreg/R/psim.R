psim <-
function (q, mu, sig) {
	ll<-length(q)
	pp<-rep(0,ll)
	dsimp <- function(x){1/sqrt(2*pi*sig^2*(x*(1-x))^3)*exp(-1/2/sig^2*(x-mu)^2/(x*(1-x)*mu^2*(1-mu)^2))}
	for (i in 1:ll) {
		if (sig < 0.001 | (1-mu)*sig < 0.01) {pp[i] <- psim.norm(q[i], mu, sig)}
		else {
			tem<-integrate(Vectorize(dsimp), lower=0, upper=q[i])
			pp[i]<-tem$value
		}
	}
	return(pp)
	}
