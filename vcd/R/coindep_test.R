coindep_test <- function(x, margin = NULL, n = 1000,
  indepfun = function(x) max(abs(x)), aggfun = max,
  alternative = c("greater", "less"),
  pearson = TRUE)
{
  DNAME <- deparse(substitute(x))
  alternative <- match.arg(alternative)

  if(is.null(margin)) {
    rs <- rowSums(x)
    cs <- colSums(x)
    expctd <- rs %o% cs / sum(rs)
    Pearson <- function(x) (x - expctd)/sqrt(expctd)      
    resids <- Pearson(x)
    
    ff <- if(is.null(aggfun)) {
      if(pearson) function(x) aggfun(indepfun(Pearson(x)))
        else function(x) aggfun(indepfun(x))
    } else {
      if(pearson) function(x) indepfun(Pearson(x))
        else function(x) indepfun(x)    
    }
       
    if(length(dim(x)) > 2) stop("currently only implemented for (conditional) 2d tables")
    dist <- sapply(r2dtable(n, rowSums(x), colSums(x)), ff)
    STATISTIC <- ff(x)   
  } else {
    ff <- if(pearson) function(x) indepfun(Pearson(x))
      else function(x) indepfun(x)

    cox <- co_table(x, margin)
    nc <- length(cox)
    if(length(dim(cox[[1]])) > 2) stop("currently only implemented for conditional 2d tables")
  
    dist <- matrix(rep(0, n * nc), ncol = nc)
    for(i in 1:nc) {
      coxi <- cox[[i]]    
      cs <- colSums(coxi)
      rs <- rowSums(coxi)
      expctd <- rs %o% cs / sum(rs)
      Pearson <- function(x) (x - expctd)/sqrt(expctd)      
      if(any(c(cs, rs) < 1)) warning("structural zeros") ## FIXME
      dist[, i] <- sapply(r2dtable(n, rs, cs), ff)
    }
    dist <- apply(dist, 1, aggfun)

    Pearson <- function(x) {
      expctd <- rowSums(x) %o% colSums(x) / sum(x)
      return((x - expctd)/sqrt(expctd))
    }
    STATISTIC <- aggfun(sapply(cox, ff))

    ## just for returning nicely formatted fitted values
    ## and residuals: fit once more with loglm()
    vars <- names(dimnames(x))
    condvars <- if(is.numeric(margin)) vars[margin] else margin
    indvars <- vars[!(vars %in% condvars)]
    coind.form <- as.formula(paste("~ (", paste(indvars, collapse = " + "),
				   ") * ",
				   paste(condvars, collapse = " * "),
				   sep = ""))
    fm <- loglm(coind.form, data = x, fitted = TRUE)
    expctd <- fitted(fm)
    resids <- residuals(fm, type = "pearson")
  }

  pdist <- function(x) sapply(x, function(y) mean(dist <= y))
  qdist <- function(p) quantile(dist, p)
  PVAL <- switch(alternative,
                 greater = mean(dist >= STATISTIC),
                 less = mean(dist <= STATISTIC))
  METHOD <- "Permutation test for conditional independence"
  names(STATISTIC) <- "f(x)"

  rval <- list(statistic = STATISTIC,
               p.value = PVAL,
	       method = METHOD,
	       data.name = DNAME,
	       observed = x,
	       expected = expctd,
	       residuals = resids,
	       margin = margin,	       
	       dist = dist,
	       qdist = qdist,
	       pdist = pdist)
  class(rval) <- c("coindep_test", "htest")
  return(rval)
}

fitted.coindep_test <- function(object, ...)
  object$expected

## plot.coindep_test
## mosaic.coindep_test
## assoc.coindep_test
## difficult, depends on functionals...
