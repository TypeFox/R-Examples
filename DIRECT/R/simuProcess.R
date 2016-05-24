simuProcess <-
function (start, beta, mu, sigma, dt, T, model=c("OU", "BM"), plot=FALSE, ...)
{
	model = match.arg (model)
	
	time = seq (0, T, dt)
	n = T/dt
	normran = rnorm (n, 0, 1)
	
	if (model=="OU")
	{
		a = exp (-beta * dt)
		b = mu * (1-a)
		c = sigma * sqrt ((1-a^2) / 2 / beta) * normran
	}
	else
	{
		a = 1
		b = 0
		c = sigma * dt * normran
	}
	
	path = rep (0, n+1)
	path[1] = start
	for (i in 2:(n+1))
		path[i] = a*path[i-1] + b + c[i-1]
	if (plot==TRUE)
		plot (time, path, type="l", xlab="time", ylab="path", ...)
	
	return (cbind (time, path))
}

