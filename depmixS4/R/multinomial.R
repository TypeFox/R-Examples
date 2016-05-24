# 
# MAKE FAMILY OBJECT FOR MULTINOMIAL RESPONSES WITH LOGISTIC LINK
# 

multinomial <-
function(link="mlogit",base=1) {
	# adapted from gaussian()
	linktemp <- substitute(link)
	if (!is.character(linktemp)) {
		linktemp <- deparse(linktemp)
		if (linktemp == "link") {
			warning("use of multinomial(link=link) is deprecated\n",
				domain = NA)
			linktemp <- eval(link)
			if (!is.character(linktemp) || length(linktemp) !=1)
			stop("'link' is invalid", domain = NA)
		}
	}
	okLinks <- c("mlogit")
	if (linktemp %in% okLinks) {
		if(linktemp == "mlogit") stats <- mlogit() else stats <- make.link(linktemp)
	} else {
		if(is.character(link)) {
			stats <- make.link(link)
			linktemp <- link
		} else {
			if(inherits(link, "link-glm")) {
				stats <- link
				if (!is.null(stats$name))
				linktemp <- stats$name
			} else {
				stop(gettextf("link \"%s\" not available for multinomial family; available links are %s",
						linktemp, paste(sQuote(okLinks), collapse = ", ")),
					domain = NA)
			}
		}
	}
	variance <- function(mu) {
		n <- length(mu)
		v <- diag(n)*outer(mu,1-mu) - (1-diag(n))*outer(mu,-mu)
	}
	validmu <- function(mu) {
		all(mu > 0) && all(mu < 1)
	}
	
	dev.resids <- function(y,mu,wt) {	
	}
	initialize <- expression()
	structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun,
			linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
			mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu, valideta = stats$valideta, base=base),
		class = "family")
}

