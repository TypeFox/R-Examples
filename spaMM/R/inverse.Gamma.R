## InverseGammaDistribution[1 + \[Nu], \[Nu] \[Mu]] (NOT linear exp family member)
## a family pseudo-definition
## see standard families in
## https://svn.r-project.org/R/trunk/src/library/stats/R/family.R
## https://svn.r-project.org/R/trunk/src/library/stats/src/family.c
## note that hglm defines `inverse.gamma`, with an "inverse" link. Not clear how they use $linkinv
`inverse.Gamma` <-function (link = "-1/mu") { ## v=-1/u, p.181, the sign of v must matter...
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    okLinks <- c("-1/mu","log")
    if (linktemp %in% okLinks)
        if (linktemp == "-1/mu") { ## case not handled by make.link but obvious modif of "inverse"
             stats<-list()
             stats$linkfun <- function(mu) -1/mu ##  v=-1/u, p.181
             stats$linkinv <- function(eta) -1/eta
             stats$mu.eta <- function(eta) 1/(eta^2)
             stats$valideta <- function(eta) all(eta!=0)
        } else stats <- make.link(linktemp)
    else {
       stop(gettextf("link \"%s\" not available for inverse-gamma family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
    }
    structure(list(family = "inverse.Gamma",
		   link = linktemp,
		   linkfun = stats$linkfun,
		   linkinv = stats$linkinv,
		   variance = function(mu) mu^2,
           ## a guessed fn of mu, but we only need it with mu=1 => LeeN01 Table 1 
		   dev.resids = function(y,mu,wt) {2 * wt * (log(ifelse(y == 0, 1, y/mu)) + (mu - y)/y)}, 
		   aic = NULL,
		   mu.eta = stats$mu.eta,
		   initialize = NULL,
		   validmu = NULL,
		   valideta = stats$valideta,
                   simulate = NULL),
	class = "family")
} 
