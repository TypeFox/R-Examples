#Copyright (c) 2015 Abby Rothway

#You're welcome to use this program or modified versions of it
#However, please cite me if you do

png = function (n)
{	n = as.integer (n) [1]
	if (n > 0)
	{	v = integer (n)
		v [1] = 2
		if (n > 1)
		{	m = 3
			for (i in 2:n)
			{	while (! is_prime (m) )
					m = m + 2
				v [i] = m
				m = m + 2
			}
		}
		v
	}
	else
		stop ("png requires positive integer")
}

is_prime = function (m)
{	m = as.integer (m)
	if (all (m > 0) )
	{	n = length (m)
		v = logical (n)
		for (i in 1:n)
			v [i] = m [i] > 1 && is.na (pnf (m [i]) )
		v
	}
	else
		stop ("is_prime requires positive integers")
}

pnf = function (m)
{	m = as.integer (m) [1]
	if (m > 0)
	{	if (m < 4 || m == 5 || m == 7)
			NA
		else if (m %% 2 == 0)
			2
		else
		{	f = NA
			x = as.integer (ceiling (ceiling (sqrt (m) ) / 2) )
			r = (1:x) * 2L + 1L
			i = 1
			n = length (r)
			while (is.na (f) && i <= n)
			{	if (m %% r [i] == 0)
					f = r [i]
				i = i + 1
			}
			f
		}
	}
	else
		stop ("pnf requires positive integer")
}

factorize = function (m)
{	m = as.integer (m) [1]
	if (m > 1 && ! is_prime (m) )
	{	i = 1
		t = pnf (m)
		r = integer ()
		while (! is.na (t) )
		{	r [i] = t
			i = i + 1
			m = m %/% t
			t = pnf (m)
		}
		r [i] = m
		f = cbind (unique (r), 0)
		for (i in 1:nrow (f) )
			f [i, 2] = sum (f [i, 1] == r)
		f
	}
	else
		stop ("can not factorize")
}
