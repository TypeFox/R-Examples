rjordan <-
function(letter, n, ...)
{
	switch(letter,
		a = rt(n, 3)/sqrt(3),
		b = rmix.dexp(n)/sqrt(2),
		c = (runif(n) - 0.5) * sqrt(12),
		d = rt(n, 5)/sqrt(5/3),
		e =  - logb(runif(n)) - 1,
		f = rmix.dexp(n, c(-3, 3))/sqrt(3^2 + 2),
		g = r.gaussmix(n, c(-2.5, 2.5)),
		h = r.gaussmix(n, c(-1.2, 1.2)),
		i = r.gaussmix(n, c(-1, 1)),
		j = r.gaussmix(n, c(-2.5, 2.5), c(0.75, 0.25)),
		k = r.gaussmix(n, c(-1.7, 1.7), c(0.75, 0.25)),
		l = r.gaussmix(n, c(-1.2, 1.2), c(0.75, 0.25)),
		m = r.gaussmix(n, c(-6, -2, 2, 6), c(0.15,
			0.35, 0.35, 
			0.15)),
		n = r.gaussmix(n, c(-4, -1, 1, 4), c(0.15,
			0.35, 0.35, 
			0.15)),
		o = r.gaussmix(n, c(-3, -0.8, 
			0.8, 3), c(0.2, 
			0.3, 0.3, 
			0.2)),
		p = r.gaussmix(n, c(-6, -2, 1, 5), c(0.2,
			0.2, 0.45,
			0.15)),
		q = r.gaussmix(n, c(-4, -1, 1, 4), c(0.1,
			0.35, 0.4, 
			0.15)),
		r = r.gaussmix(n, c(-3, -1, 0.8, 3.5), c(
			0.1, 0.35, 
			0.4, 0.15)))
}

