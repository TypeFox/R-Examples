`Beta` <-
function(link="logit"){

 linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)          
        if (linktemp == "link") {
            warning("use of beta(link=link) is deprecated\n",
                domain = NA)
            linktemp <- eval(link)
            if (!is.character(linktemp) || length(linktemp) !=
                1)
                stop("'link' is invalid", domain = NA)
        }
    }
    okLinks <- c("logit","probit","cloglog","log", "identity")
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
            stop(gettextf("link \"%s\" not available for beta family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }

  variance<-function(arg1){
    return(arg1*(1-arg1))
  }
  dev.resids<-function(y, mu, wt) {
    dev <- 2 * wt * (y*log( y/mu) +(1-y)*log((1-y)/(1-mu)))
    return(pmax(dev,0)) ## resids(1/2,1/2,1) might be slightly negative...
  }

structure(list(family = "Beta", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, mu.eta = stats$mu.eta),
        class = "family")
}

