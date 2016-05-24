.dsgt = function(x, mu, sigma, lambda, p, q, mean.cent, var.adj, log) {
	if(var.adj) sigma = sigma / (q^(1/p)*sqrt((3*lambda^2+1)*(beta(3/p, q-2/p)/beta(1/p,q))-4*lambda^2*(beta(2/p, q-1/p)/beta(1/p, q))^2))
	if(mean.cent) x = x + (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	if(!log) return(p/(2*sigma*q^(1/p)*beta(1/p,q)*(1+abs(x-mu)^p/(q*sigma^p*(1+lambda*sgn(x-mu))^p))^(q+1/p)))
	return(log(p)-log(2)-log(sigma)-log(q)/p-lbeta(1/p,q)-(1/p+q)*log(1+abs(x-mu)^p/(q*sigma^p*(1+lambda*sgn(x-mu))^p)))
}

.psgt = function(quant, mu, sigma, lambda, p, q, mean.cent, var.adj, lower.tail, log.p) {
  if(var.adj) sigma = sigma / (q^(1/p)*sqrt((3*lambda^2+1)*(beta(3/p, q-2/p)/beta(1/p,q))-4*lambda^2*(beta(2/p, q-1/p)/beta(1/p, q))^2))
  if(mean.cent) quant = quant + (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	quant = quant - mu
	flip = quant > 0
	lambda[flip] = -lambda[flip]
	quant[flip] = -quant[flip]
	out = (1-lambda)/2+(lambda-1)/2*stats::pbeta(1/(1+q*(sigma*(1-lambda)/(-quant))^p),1/p,q)
	out[flip] = 1 - out[flip]
	if(!lower.tail) out = 1 - out
	if(log.p) out = log(out)
	return(out)
}

.qsgt = function(prob, mu, sigma, lambda, p, q, mean.cent, var.adj, lower.tail, log.p) {
  if(var.adj) sigma = sigma / (q^(1/p)*sqrt((3*lambda^2+1)*(beta(3/p, q-2/p)/beta(1/p,q))-4*lambda^2*(beta(2/p, q-1/p)/beta(1/p, q))^2))
  if(log.p) prob = exp(prob)
	if(!lower.tail)	prob = 1 - prob
	flip = prob > (1-lambda)/2
	prob[flip] = 1 - prob[flip]
	lam = lambda
	lam[flip] = -lam[flip]
	out = sigma*(lam-1)*(1/(q*stats::qbeta(1-2*prob/(1-lam),1/p,q))-1/q)^(-1/p)
	out[flip] = -out[flip]
	out = out + mu
	if(mean.cent) out = out - (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	return(out)
}
