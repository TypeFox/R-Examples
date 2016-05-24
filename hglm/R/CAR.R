`CAR` <-
function(D, link = "identity", link.rand.disp = "inverse"){
linktemp <- substitute(link)
linktemp.rand.disp <- substitute(link.rand.disp)
if (!is.character(linktemp)) 
	linktemp <- deparse(linktemp)
if (!is.character(linktemp.rand.disp)) 
	linktemp.rand.disp <- deparse(linktemp.rand.disp)
okLinks <- c("inverse", "log", "identity")
okLinks.rand.disp <- c("inverse", "log")
if (linktemp %in% okLinks) 
	stats <- make.link(linktemp)
else if (is.character(link)) {
	stats <- make.link(link)
	linktemp <- link
}
else {
	if (inherits(link, "link-glm")) {
		stats <- link
		if (!is.null(stats$name)) 
			linktemp <- stats$name
	}
	else {
		stop(gettextf("link \"%s\" not available for CAR family; available links are %s", 
						linktemp, paste(sQuote(okLinks), collapse = ", ")), 
				domain = NA)
	}
}
if (linktemp.rand.disp %in% okLinks.rand.disp) 
	stats.rand.disp <- make.link(linktemp.rand.disp)
else if (is.character(link.rand.disp)) {
	stats.rand.disp <- make.link(link.rand.disp)
	linktemp.rand.disp <- link.rand.disp
}
else {
	if (inherits(link.rand.disp, "link-glm")) {
		stats.rand.disp <- link.rand.disp
		if (!is.null(stats.rand.disp$name)) 
			linktemp.rand.disp <- stats.rand.disp$name
	}
	else {
		stop(gettextf("link.rand.disp \"%s\" not available for CAR family; available links are %s", 
						linktemp.rand.disp, paste(sQuote(okLinks.rand.disp), collapse = ", ")), 
				domain = NA)
	}
}
decomp <- eigen(D)
structure(list(family = "CAR", 
				link = linktemp, link.rand.disp = linktemp.rand.disp, 
				linkfun = stats$linkfun, linkfun.rand.disp = stats.rand.disp$linkfun, 
				linkinv = stats$linkinv, linkinv.rand.disp = stats.rand.disp$linkinv, 
				Dvec = decomp$vectors, Deigen = decomp$values,
				variance = function(mu) rep.int(1, length(mu)), 
				dev.resids = function(y, mu, wt) wt * ((y - mu)^2), 
				aic = function(y, n, mu, wt, dev) {
					nobs <- length(y)
					nobs * (log(dev/nobs * 2 * pi) + 1) + 2 - sum(log(wt))
				}, 
				mu.eta = stats$mu.eta, 
				initialize = expression({
							n <- rep.int(1, nobs)
							if (is.null(etastart) && is.null(start) && is.null(mustart) && 
									((family$link == "inverse" && any(y == 0)) || 
										(family$link == "log" && any(y <= 0)))) stop("cannot find valid starting values: please specify some")
							mustart <- y
						}), validmu = function(mu) TRUE, valideta = stats$valideta), 
		class = "family")
}

