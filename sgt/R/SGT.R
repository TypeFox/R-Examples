dsgt = function(x, mu = 0, sigma = 1, lambda = 0, p = 2, q = Inf, mean.cent = TRUE, var.adj = TRUE, log = FALSE) {
	n = max(length(x), length(mu), length(sigma), length(lambda), length(p), length(q))
	x = rep(x, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	out = rep(NaN, length.out=n)
	if(!is.logical(var.adj[1L])) {
	  if(is.finite(var.adj[1L]) & var.adj[1L] > 0) { sigma = sigma * var.adj
	  } else {warning("var.adj argument is invalid and assumed FALSE")}
	  var.adj = FALSE
	}
	goodarg = is.finite(mu) & is.finite(sigma) & sigma > 0 & lambda < 1 & lambda > -1 & p > 0 & q > 0
	goodmean = !mean.cent[1L] | p*q > 1
	goodvar = !var.adj[1L] | p*q > 2
	goodall = goodarg & goodmean & goodvar
	if(sum(goodarg) != n) warning("NaNs produced")
	if(sum(goodmean) != n) warning("mean.cent is TRUE but p*q <= 1; NaNs produced")
	if(sum(goodmean) != n) warning("var.adj is TRUE but p*q <= 2; NaNs produced")
  	goodsgt = goodall & is.finite(q) & is.finite(p)
  	out[goodsgt] = .dsgt(x[goodsgt], mu[goodsgt], sigma[goodsgt], lambda[goodsgt], p[goodsgt], q[goodsgt], mean.cent[1L], var.adj[1L], log[1L])
  	goodsged = goodall & is.finite(p) & q == Inf
  	out[goodsged] = .dsged(x[goodsged], mu[goodsged], sigma[goodsged], lambda[goodsged], p[goodsged], mean.cent[1L], var.adj[1L], log[1L])
	goodunif = goodall & p == Inf
  	out[goodunif] = .dunif(x[goodunif], mu[goodunif], sigma[goodunif], var.adj[1L], log[1L])
  	return(out)
}

psgt = function(quant, mu = 0, sigma = 1, lambda = 0, p = 2, q = Inf, mean.cent = TRUE, var.adj = TRUE, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(quant), length(mu), length(sigma), length(lambda), length(p), length(q))
	quant = rep(quant, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	out = rep(NaN, length.out=n)
	if(!is.logical(var.adj[1L])) {
	  if(is.finite(var.adj[1L]) & var.adj[1L] > 0) { sigma = sigma * var.adj
	  } else {warning("var.adj argument is invalid and assumed FALSE")}
	  var.adj = FALSE
	}
	goodarg = is.finite(mu) & is.finite(sigma) & sigma > 0 & lambda < 1 & lambda > -1 & p > 0 & q > 0
	goodmean = !mean.cent[1L] | p*q > 1
	goodvar = !var.adj[1L] | p*q > 2
	goodall = goodarg & goodmean & goodvar
	if(sum(goodarg) != n) warning("NaNs produced")
	if(sum(goodmean) != n) warning("mean.cent is TRUE but p*q <= 1; NaNs produced")
	if(sum(goodmean) != n) warning("var.adj is TRUE but p*q <= 2; NaNs produced")
  	goodsgt = goodall & is.finite(q) & is.finite(p) & is.finite(quant)
  	out[goodsgt] = .psgt(quant[goodsgt], mu[goodsgt], sigma[goodsgt], lambda[goodsgt], p[goodsgt], q[goodsgt], mean.cent[1L], var.adj[1L], lower.tail[1L], log.p[1L])
  	goodsged = goodall & is.finite(p) & q == Inf & is.finite(quant)
  	out[goodsged] = .psged(quant[goodsged], mu[goodsged], sigma[goodsged], lambda[goodsged], p[goodsged], mean.cent[1L], var.adj[1L], lower.tail[1L], log.p[1L])
	goodunif = goodall & p == Inf & is.finite(quant)
  	out[goodunif] = .punif(quant[goodunif], mu[goodunif], sigma[goodunif], var.adj[1L], lower.tail[1L], log.p[1L])
  	out[goodall & quant == -Inf] = 0
  	out[goodall & quant == Inf] = 1
  	return(out)
}

qsgt = function(prob, mu = 0, sigma = 1, lambda = 0, p = 2, q = Inf, mean.cent = TRUE, var.adj = TRUE, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(prob), length(mu), length(sigma), length(lambda), length(p), length(q))
	prob = rep(prob, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	out = rep(NaN, length.out=n)
	if(!is.logical(var.adj[1L])) {
	  if(is.finite(var.adj[1L]) & var.adj[1L] > 0) { sigma = sigma * var.adj
	  } else {warning("var.adj argument is invalid and assumed FALSE")}
	  var.adj = FALSE
	}
	goodarg = is.finite(mu) & is.finite(sigma) & sigma > 0 & lambda < 1 & lambda > -1 & p > 0 & q > 0
	goodmean = !mean.cent[1L] | p*q > 1
	goodvar = !var.adj[1L] | p*q > 2
	goodall = goodarg & goodmean & goodvar
	if(sum(goodarg) != n) warning("NaNs produced")
	if(sum(goodmean) != n) warning("mean.cent is TRUE but p*q <= 1; NaNs produced")
	if(sum(goodmean) != n) warning("var.adj is TRUE but p*q <= 2; NaNs produced")
  	goodsgt = goodall & is.finite(q) & is.finite(p)
  	out[goodsgt] = .qsgt(prob[goodsgt], mu[goodsgt], sigma[goodsgt], lambda[goodsgt], p[goodsgt], q[goodsgt], mean.cent[1L], var.adj[1L], lower.tail[1L], log.p[1L])
  	goodsged = goodall & is.finite(p) & q == Inf
  	out[goodsged] = .qsged(prob[goodsged], mu[goodsged], sigma[goodsged], lambda[goodsged], p[goodsged], mean.cent[1L], var.adj[1L], lower.tail[1L], log.p[1L])
	goodunif = goodall & p == Inf
  	out[goodunif] = .qunif(prob[goodunif], mu[goodunif], sigma[goodunif], var.adj[1L], lower.tail[1L], log.p[1L])
  	return(out)
}

rsgt = function(n, mu = 0, sigma = 1, lambda = 0, p = 2, q = Inf, mean.cent = TRUE, var.adj = TRUE) {
	if(length(n) > 1) n = length(n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	return(qsgt(stats::runif(n), mu, sigma, lambda, p, q, mean.cent, var.adj))
}
