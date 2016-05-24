vonMises <- function (link = "tan") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("tan", "log", "probit", "identity")
  if (linktemp %in% okLinks)
    stats <- make.circular.link(linktemp)
  else if (is.character(link)) {
    stats <- make.circular.link(link)
    linktemp <- link
  } else {
  ## what else shall we allow?  At least objects of class link-glm.
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for gaussian family; available links are %s',
	linktemp, paste(sQuote(okLinks), collapse =", ")),
	domain = NA)
    }
  }
  structure(list(family = "vonMises",
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    variance = function(mu) rep.int(1, length(mu)),
    dev.resids = function(y, mu, mulinear, kappa, wt) {
      if (kappa < 100000)
        llik <- 2*sum(wt)*(log(2*pi)+log(besselI(kappa,nu=0,expon.scaled=TRUE))+kappa) -2*sum(wt * kappa * cos(y-mu-mulinear))
      else
        llik <- ifelse((y-mu-mulinear)==0, -Inf, Inf)
      return(llik)
    },
    aic = function(dev,rank){
      dev + 2*rank
    },
    mu.eta = stats$mu.eta,
    #initialize = expression({
    #  n <- rep.int(1, nobs)
    #  if(is.null(etastart) && is.null(start) &&
    #    is.null(mustart) && (family$link == "log" && any(y <= 0)))
    #    stop("cannot find valid starting values: please specify some")
    #    mustart <- y - circular:::MeanCircularRad(y)}),
    validmu = function(mu) TRUE,
    valideta = stats$valideta
  ),
  class = "family")
}
