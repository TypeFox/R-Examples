SIGCONST = 1.78732427093276085059404775102

dihs = function(x, mu = 0, sigma = SIGCONST, lambda = 0, k = 1, log = FALSE) {
	n = max(length(x), length(mu), length(sigma), length(lambda), length(k))
	x = rep(x, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	k = rep(k, length.out=n)
	b = 2*sigma / sqrt((exp(2*lambda+k^(-2))+exp(-2*lambda+k^(-2))+2)*(exp(k^(-2))-1))
    a = mu - 0.5*b*(exp(lambda)-exp(-lambda))*exp(0.5*k^(-2))
    if(!log[1L]) return(k*exp((-k^2/2)*(log((x-a)/b+sqrt(((x-a)/b)^2+1))-lambda)^2)/sqrt(2*pi*((a-x)^2+b^2)))
    return(log(k)+(-k^2/2)*(log((x-a)/b+sqrt(((x-a)/b)^2+1))-lambda)^2-0.5*log(2*pi*((a-x)^2+b^2)))
}

pihs = function(q, mu = 0, sigma = SIGCONST, lambda = 0, k = 1, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(q), length(mu), length(sigma), length(lambda), length(k))
	q = rep(q, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	k = rep(k, length.out=n)
	b = 2*sigma / sqrt((exp(2*lambda+k^(-2))+exp(-2*lambda+k^(-2))+2)*(exp(k^(-2))-1))
    a = mu - 0.5*b*(exp(lambda)-exp(-lambda))*exp(0.5*k^(-2))
    out = pnorm(k*log((q-a)/b+sqrt(((q-a)/b)^2+1))-lambda*k)
	if(!lower.tail[1L]) out = 1 - out
	if(log.p[1L]) out = log(out)
	return(out)
}

qihs = function(p, mu = 0, sigma = SIGCONST, lambda = 0, k = 1, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(p), length(mu), length(sigma), length(lambda), length(k))
	p = rep(p, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	k = rep(k, length.out=n)
	if(log.p[1L]) p = exp(p)
	if(!lower.tail[1L])	p = 1 - p
	b = 2*sigma / sqrt((exp(2*lambda+k^(-2))+exp(-2*lambda+k^(-2))+2)*(exp(k^(-2))-1))
    a = mu - 0.5*b*(exp(lambda)-exp(-lambda))*exp(0.5*k^(-2))
    a + b*sinh(lambda+qnorm(p)/k)
}

rihs = function(n, mu = 0, sigma = SIGCONST, lambda = 0, k = 1) {
	if(length(n) > 1) n = length(n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	k = rep(k, length.out=n)
    b = 2*sigma / sqrt((exp(2*lambda+k^(-2))+exp(-2*lambda+k^(-2))+2)*(exp(k^(-2))-1))
    a = mu - 0.5*b*(exp(lambda)-exp(-lambda))*exp(0.5*k^(-2))
    a + b*sinh(lambda+rnorm(n)/k)
}