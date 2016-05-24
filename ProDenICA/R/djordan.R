djordan <-
function(letter, x, ...)
{
	switch(letter,
		a = dt(x * sqrt(3), 3),
		b = dmix.dexp(x * sqrt(2)),
		c = dunif(x/sqrt(12) + 0.5),
		d = dt(x * sqrt(5/3), 5),
		e = exp( - (x + 1)) * (x > -1),
		f = dmix.dexp(x * sqrt(3^2 + 2), c(-3, 3)),
		g = d.gaussmix(x, c(-2.5, 2.5)),
		h = d.gaussmix(x, c(-1.2, 1.2)),
		i = d.gaussmix(x, c(-1, 1)),
		j = d.gaussmix(x, c(-2.5, 2.5), c(0.75, 0.25)),
		k = d.gaussmix(x, c(-1.7, 1.7), c(0.75, 0.25)),
		l = d.gaussmix(x, c(-1.2, 1.2), c(0.75, 0.25)),
		m = d.gaussmix(x, c(-6, -2, 2, 6), c(0.15,
			0.35, 0.35, 
			0.15)),
		n = d.gaussmix(x, c(-4, -1, 1, 4), c(0.15,
			0.35, 0.35, 
			0.15)),
		o = d.gaussmix(x, c(-3, -0.8, 
			0.8, 3), c(0.2, 
			0.3, 0.3, 
			0.2)),
		p = d.gaussmix(x, c(-6, -2, 1, 5), c(0.2,
			0.2, 0.45, 
			0.15)),
		q = d.gaussmix(x, c(-4, -1, 1, 4), c(0.1,
			0.35, 0.4, 
			0.15)),
		r = d.gaussmix(x, c(-3, -1, 0.8, 3.5), c(
			0.1, 0.35, 
			0.4, 0.15)))
}

