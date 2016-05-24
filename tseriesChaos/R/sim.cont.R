#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
sim.cont <- function(syst, start.time, end.time, dt, start.x, parms=NULL, obs.fun=function(x) x[1] ) {
	times <- seq(start.time, end.time, by=dt)
	series <- lsoda(start.x, times, func=syst, parms=parms)[,-1]
	series <- apply(series, 1, obs.fun)
	series <- ts(series, start=start.time, end=end.time, deltat=dt)
}
