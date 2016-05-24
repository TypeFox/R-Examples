sgn = function(x) {
	x[x >= 0] = 1
	x[x < 0] = -1
	return(x)
}

.dsged = function(x, mu, sigma, lambda, p, mean.cent, var.adj, log) {
  if(var.adj) sigma = sigma / sqrt((pi*(1+3*lambda^2)*gamma(3/p)-16^(1/p)*lambda^2*(gamma(1/2+1/p))^2*gamma(1/p))/(pi*gamma(1/p)))
  if(mean.cent) x = x + (2^(2/p)*sigma*lambda*gamma(1/2+1/p))/sqrt(pi)
	if(!log) return(p/(2*sigma*gamma(1/p)*exp((abs(x-mu)/(sigma*(1+lambda*sgn(x-mu))))^p)))
	return(log(p)-log(2)-log(sigma)-lgamma(1/p)-(abs(x-mu)/(sigma*(1+lambda*sgn(x-mu))))^p)
}

.psged = function(quant, mu, sigma, lambda, p, mean.cent, var.adj, lower.tail, log.p) {
  if(var.adj) sigma = sigma / sqrt((pi*(1+3*lambda^2)*gamma(3/p)-16^(1/p)*lambda^2*(gamma(1/2+1/p))^2*gamma(1/p))/(pi*gamma(1/p)))
  if(mean.cent) quant = quant + (2^(2/p)*sigma*lambda*gamma(1/2+1/p))/sqrt(pi)
	quant = quant - mu
	flip = quant < 0
	lambda[flip] = -lambda[flip]
	quant[flip] = -quant[flip]
	out = (1-lambda)/2+(1+lambda)/2*stats::pgamma((quant/(sigma*(1+lambda)))^p, 1/p)
	out[flip] = 1 - out[flip]
	if(!lower.tail) out = 1 - out
	if(log.p) out = log(out)
	return(out)
}

.qsged = function(prob, mu, sigma, lambda, p, mean.cent, var.adj, lower.tail, log.p) {
  if(var.adj) sigma = sigma / sqrt((pi*(1+3*lambda^2)*gamma(3/p)-16^(1/p)*lambda^2*(gamma(1/2+1/p))^2*gamma(1/p))/(pi*gamma(1/p)))
  if(log.p) prob = exp(prob)
	if(!lower.tail)	prob = 1 - prob
	flip = prob < (1-lambda)/2
	prob[flip] = 1 - prob[flip]
	lam = lambda
	lam[flip] = -lam[flip]
	out = sigma*(1+lam)*(stats::qgamma(2*prob/(1+lam)+(lam-1)/(lam+1), 1/p))^(1/p)
	out[flip] = -out[flip]
	out = out + mu
	if(mean.cent) out = out - (2^(2/p)*sigma*lambda*gamma(1/2+1/p))/sqrt(pi)
	return(out)
}
