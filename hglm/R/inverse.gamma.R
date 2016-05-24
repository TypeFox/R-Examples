`inverse.gamma` <-
function(link="inverse"){
 linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") {
            warning("use of inverse.gamma(link=link) is deprecated\n",
                domain = NA)
            linktemp <- eval(link)
            if (!is.character(linktemp) || length(linktemp) !=
                1)
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
            stop(gettextf("link \"%s\" not available for inverse-gamma family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }

  variance<-function(arg1){
    return(arg1^2)
  }
  dev.resids<-function(y, mu, wt) {
    dev <- -2 * wt * (log( y/mu) - (y - mu)/mu)
    return(dev)
  }

structure(list(family = "inverse.gamma", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, mu.eta = stats$mu.eta),
        class = "family")
}

