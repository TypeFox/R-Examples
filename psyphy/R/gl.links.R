logit.2asym <- function(g, lam) {
	if ((g < 0 ) || (g > 1))
		stop("g must in (0, 1)")
	if ((lam < 0) || (lam > 1))
		stop("lam outside (0, 1)")
	linkfun <- function(mu) {
		mu <- pmin(mu, 1 - (lam + .Machine$double.eps))
		mu <- pmax(mu, g + .Machine$double.eps)
		qlogis((mu - g)/(1 - g - lam))
		}
	linkinv <- function(eta) {
		g + (1 - g - lam) * binomial()$linkinv(eta)
#		 .Call("logit_linkinv", eta, PACKAGE = "stats")
		}
	mu.eta <- function(eta) {
		(1 - g - lam) * binomial()$mu.eta(eta)
#		.Call("logit_mu_eta", eta, PACKAGE = "stats")
		}
	valideta <- function(eta) TRUE
	link <- paste("logit.2asym(", g, ", ", lam, ")", sep = "")
	structure(list(linkfun = linkfun, linkinv = linkinv, 
	mu.eta = mu.eta, valideta = valideta, name = link), 
	class = "link-glm")
}

probit.2asym <- function(g, lam) {
	if ((g < 0 ) || (g > 1))
		stop("g must in (0, 1)")
	if ((lam < 0) || (lam > 1))
		stop("lam outside (0, 1)")
	linkfun <- function(mu) {
		mu <- pmin(mu, 1 - (lam + .Machine$double.eps))
		mu <- pmax(mu, g + .Machine$double.eps)
		qnorm((mu - g)/(1 - g - lam))
		}
	linkinv <- function(eta) {
		g + (1 - g - lam) * 
		 pnorm(eta)
		}
	mu.eta <- function(eta) {
		(1 - g - lam) * dnorm(eta)		}
	valideta <- function(eta) TRUE
	link <- paste("probit.2asym(", g, ", ", lam, ")", sep = "")
	structure(list(linkfun = linkfun, linkinv = linkinv, 
	mu.eta = mu.eta, valideta = valideta, name = link), 
	class = "link-glm")
}

cauchit.2asym <- function(g, lam) {
	if ((g < 0 ) || (g > 1))
		stop("g must in (0, 1)")
	if ((lam < 0) || (lam > 1))
		stop("lam outside (0, 1)")
	linkfun <- function(mu) {
		mu <- pmin(mu, 1 - (lam + .Machine$double.eps))
		mu <- pmax(mu, g + .Machine$double.eps)
		qcauchy((mu - g)/(1 - g - lam))
		}
	linkinv <- function(eta) {
		g + (1 - g - lam) * 
		 pcauchy(eta)
		}
	mu.eta <- function(eta) {
		(1 - g - lam) * dcauchy(eta)		}
	valideta <- function(eta) TRUE
	link <- paste("cauchit.2asym(", g, ", ", lam, ")", sep = "")
	structure(list(linkfun = linkfun, linkinv = linkinv, 
	mu.eta = mu.eta, valideta = valideta, name = link), 
	class = "link-glm")
}

cloglog.2asym <- function(g, lam) {
	if ((g < 0 ) || (g > 1))
		stop("g must in (0, 1)")
	if ((lam < 0) || (lam > 1))
		stop("lam outside (0, 1)")
	linkfun <- function(mu) {
        	mu <- pmax(pmin(mu, 1 - (lam + .Machine$double.eps)), 
        		g + .Machine$double.eps)
        	log(-log((mu - g)/(1 - g - lam)))
    	}
	linkinv <- function(eta) {
       	 tmp <- g + (1 - g - lam) * (-expm1(-exp(eta)))
         pmax(pmin(tmp, 1 - (lam + .Machine$double.eps)), 
        	g + .Machine$double.eps)
    	}
	mu.eta <- function(eta) {
          eta <- pmin(eta, 700)
          pmax((1 - g - lam) * exp(eta) * exp(-exp(eta)), 
        	.Machine$double.eps)
	    }
	valideta <- function(eta) TRUE
	link <- paste("cloglog.2asym(", g, ", ", lam, ")", sep = "")
	structure(list(linkfun = linkfun, linkinv = linkinv, 
	mu.eta = mu.eta, valideta = valideta, name = link), 
	class = "link-glm")
}

weib.2asym <- function(...) cloglog.2asym(...)	
	
	
