quasibinomialMany <-
function (link = "logit")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") {
            warning("use of quasibinomial(link=link) is deprecated\n",
                domain = NA)
            linktemp <- eval(link)
            if (!is.character(linktemp) || length(linktemp) !=
                1)
                stop("'link' is invalid", domain = NA)
        }
    }
    okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
    if(linktemp=="varstab"){
        stats <- make.vslink("quasibinomial")
        linktemp <- stats$name    # however this hasn't an attribute, as there is no easy analogue
    } else if (linktemp %in% okLinks) {
        stats <- make.link(linktemp)
    } else if (is.character(link)) {
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
            stop(gettextf("link \"%s\" not available for quasibinomial family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(mu > 0) && all(mu < 1)

    dev.resids <- function(y, mu, wt) {
      y <- as.matrix(y)
      mu <- as.matrix(mu)
      q <- ncol(y)
      dr <- matrix(ncol=q, nrow=nrow(y))
      for(i in 1:q)
        dr[,i] <- 2 * wt * (y[,i] * log(ifelse(y[,i] ==
        0, 1, y[,i]/mu[,i])) + (1 - y[,i]) * log(ifelse(y[,i] == 1,
        1, (1 - y[,i])/(1 - mu[,i]))))
      return(dr)
    }
    aic <- function(y, n, mu, wt, dev) rep(NA, times=NCOL(y))
    
    initialize <- expression({
    
      if(is.null(colnames(y))) { berno <- TRUE
        } else if( all(substr(colnames(y), 1,4) %in% c("succ", "fail") )
           # & !is.na(as.numeric(substr(colnames(y), 2,2)))
           ){
            berno <- FALSE
        } else berno <- TRUE
        
        if ( berno ) {
            if (is.factor(y)) y <- y != levels(y)[1]
            n <- matrix(1, nrow=nobs, ncol=NCOL(y))
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            
        } else {
            struc <- substr(colnames(y),1,4)
            n <- y[, struc=="succ", drop=FALSE] +
              y[, struc=="fail", drop=FALSE]
            # if( any(diag(cov(t(n)))!=0))
            if(any(n - n[,1] != 0))
                stop("the number of trials must be the same for each variable")
            y <- ifelse(n == 0, 0, y[, struc=="succ", drop=FALSE]/n)
            weights <- n[,1] * weights    # this is a matrix --> errors possible
            mustart <- (n * y + 0.5)/(n + 1)
            
        } # else stop("for the quasibinomialMany family, y must be a matrix of 0 and 1's\n",
          #  "or a matrix with columns giving no. of successes ",
			    #  "and corresponding columns with no. of failures\n",
          #  "and the matrix' structure described by the ",
			    #  "column names (see 'binstruc')")
    })
    structure(list(family = "quasibinomial", link = linktemp,
        linkfun = stats$linkfun, linkinv = stats$linkinv, variance = variance,
        dev.resids = dev.resids, aic = aic, mu.eta = stats$mu.eta,
        initialize = initialize, validmu = validmu, valideta = stats$valideta),
        class = c("family.mvabund", "family"))
}


