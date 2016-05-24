`GAMMA` <- 
function (link = "log") { ### Modified Gamma() function ###
	linktemp <- substitute(link)
	if (!is.character(linktemp)) {
		linktemp <- deparse(linktemp)
		if (linktemp == "link") {
			warning("use of GAMMA(link=link) is deprecated\n", domain = NA)
			linktemp <- eval(link)
			if (!is.character(linktemp) || length(linktemp) != 1)
				stop("'link' is invalid", domain = NA)
		}
	}
	okLinks <- c("inverse", "log", "identity")
	if (linktemp %in% okLinks)
		stats <- make.link(linktemp)
	else if (is.character(link))
		stats <- make.link(link)
	else {
		if (inherits(link, "link-glm")) {
			stats <- link
			if (!is.null(stats$name))
				linktemp <- stats$name
		}
		else {
			stop(gettextf("link \"%s\" not available for gamma family; available links are %s",
							linktemp, paste(sQuote(okLinks), collapse = ", ")), domain = NA)
		}
	}
	variance <- function(mu) mu ### Note this variance function
	validmu <- function(mu) all(mu > 0)
	dev.resids <- function(y, mu, wt) -2 * wt * (log(mu) + (1 - mu))
	structure(list(family = "GAMMA", link = linktemp, linkfun = stats$linkfun,
					linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
					mu.eta = stats$mu.eta), class = "family")
}
