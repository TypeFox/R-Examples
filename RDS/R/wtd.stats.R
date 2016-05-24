# # From Hmisc 3.16-0 by Frank E Harrell Jr, with contributions from Charles Dupont and many others.
wtd.mean <- function(x, weights=NULL, normwt='ignored', na.rm=TRUE)
{
  if(!length(weights)) return(base::mean(x, na.rm=na.rm))
  if(na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }

  sum(weights*x)/sum(weights)
}

wtd.var <- function(x, weights=NULL, normwt=FALSE, na.rm=TRUE,
                    method = c('unbiased', 'ML'))
{
  method <- match.arg(method)
  if(!length(weights)) {
    if(na.rm) x <- x[!is.na(x)]
    return(stats::var(x))
  }

  if(na.rm) {
    s       <- !is.na(x + weights)
    x       <- x[s]
    weights <- weights[s]
  }

  if(normwt)
    weights <- weights * length(x) / sum(weights)

  if(method == 'ML')
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = "ML")$cov))

  sw   <- sum(weights)
  xbar <- sum(weights * x) / sw
  sum(weights*((x - xbar)^2)) /
    (sw - (if(normwt) sum(weights ^ 2) / sw else 1))
}


wtd.quantile <- function(x, weights=NULL, probs=c(0, .25, .5, .75, 1), 
                         normwt=FALSE, na.rm=TRUE)
{
  if(!length(weights))
    return(stats::quantile(x, probs=probs, na.rm=na.rm))

  if(any(probs < 0 | probs > 1))
    stop("Probabilities must be between 0 and 1 inclusive")

  nams <- paste(format(round(probs * 100, if(length(probs) > 1) 
                             2 - log10(diff(range(probs))) else 2)), 
                "%", sep = "")

    w <- wtd.table(x, weights, na.rm=na.rm, normwt=normwt)
    x     <- w$x
    wts   <- w$sum.of.weights
    n     <- sum(wts)
    order <- 1 + (n - 1) * probs
    low   <- pmax(floor(order), 1)
    high  <- pmin(low + 1, n)
    order <- order %% 1
    ## Find low and high order statistics
    ## These are minimum values of x such that the cum. freqs >= c(low,high)
    allq <- approx(cumsum(wts), x, xout=c(low,high), 
                   method='constant', f=1, rule=2)$y
    k <- length(probs)
    quantiles <- (1 - order)*allq[1:k] + order*allq[-(1:k)]
    names(quantiles) <- nams
    return(quantiles)
}

wtd.table <- function(x, weights=NULL, 
                      normwt=FALSE, na.rm=TRUE)
{
  if(!length(weights))
    weights <- rep(1, length(x))

  ax <- attributes(x)
  ax$names <- NULL
  
  if(is.character(x)) x <- as.factor(x)
  lev <- levels(x)
  x <- unclass(x)
  
  if(na.rm) {
    s <- !is.na(x + weights)
    x <- x[s, drop=FALSE]    ## drop is for factor class
    weights <- weights[s]
  }

  n <- length(x)
  if(normwt)
    weights <- weights * length(x) / sum(weights)

  i <- order(x)  # R does not preserve levels here
  x <- x[i]; weights <- weights[i]

  if(anyDuplicated(x)) {  ## diff(x) == 0 faster but doesn't handle Inf
    weights <- tapply(weights, x, sum)
    if(length(lev)) {
      levused <- lev[sort(unique(x))]
      if((length(weights) > length(levused)) &&
         any(is.na(weights)))
        weights <- weights[!is.na(weights)]

      if(length(weights) != length(levused))
        stop('program logic error')

      names(weights) <- levused
    }

    if(!length(names(weights)))
      stop('program logic error')

    x <- all.is.numeric(names(weights))

    names(weights) <- NULL
    return(list(x=x, sum.of.weights=weights))
  }

  xx <- x

    list(x=if(length(lev))lev[x]
           else xx, 
         sum.of.weights=weights)
}
all.is.numeric <- function(x, extras=c('.','NA'))
{
  x <- sub('[[:space:]]+$', '', x)
  x <- sub('^[[:space:]]+', '', x)
  xs <- x[match(x,c('',extras),nomatch = 0) == 0]
  isnum <- suppressWarnings(!any(is.na(as.numeric(xs))))
  if(isnum){as.numeric(x)}else{x}
}
