`var.periodic.correlation` <-
function(l, m, n, p)
{
	if((l %% p) == 0)
		(n - l/p)/(n * (n + 2))
	else (n - trunc((l - m + p)/p))/n^2
}

