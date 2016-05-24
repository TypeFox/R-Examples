"dgev" <- 
function(x, xi = 1, mu = 0, sigma = 1)
{
	tmp <- (1 + (xi * (x - mu))/sigma)
	(as.numeric(tmp > 0) * (tmp^(-1/xi - 1) * exp( - tmp^(-1/xi))))/
		sigma
}

"dgpd" <- 
function(x, xi, mu = 0, beta = 1)
{
	(beta^(-1)) * (1 + (xi * (x - mu))/beta)^((-1/xi) - 1)
}

"pgev" <- 
function(q, xi = 1, mu = 0, sigma = 1)
{
	exp( - (1  + (xi * (q - mu))/sigma)^(-1 /xi))
}

"pgpd" <- 
function(q, xi, mu = 0, beta = 1)
{
	(1  - (1  + (xi * (q - mu))/beta)^(-1 /xi))
}

"qgev" <- 
function(p, xi = 1, mu = 0, sigma = 1)
{
	mu + (sigma/xi) * (( - logb(p))^( - xi) - 1)
}

"qgpd" <- 
function(p, xi, mu = 0, beta = 1)
{
	mu + (beta/xi) * ((1 - p)^( - xi) - 1)
}

"rgev" <- 
function(n, xi = 1, mu = 0, sigma = 1)
{
	mu - (sigma * (1 - ( - logb(runif(n)))^( - xi)))/xi
}

"rgpd" <- 
function(n, xi, mu = 0, beta = 1)
{
	mu + (beta/xi) * ((1 - runif(n))^( - xi) - 1)
}
