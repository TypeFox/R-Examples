
## ########################################################
##
## Testing for condional independence in MVN-distribution
## <x>  : list(cov=, n.obs=)
## <set>: NULL, a vector or a formula
##
## ########################################################

ciTest_mvn <- function(x,set=NULL, statistic="DEV",...){
  if(any(is.na(match(c("cov","n.obs"), names(x))))){
    stop("Expecting a list with components 'cov' and 'n.obs'\n")
  }

  if (is.null(set)){
    set <- colnames(x$cov)
    x$cov <- x$cov[set,set]
  } else {
    if (inherits(set,c("formula","character"))){
      set <- unlist(rhsFormula2list(set))
      set <- colnames(x$cov)[pmatch(set, colnames(x$cov))]
      x$cov <- x$cov[set,set]
    }
  }
  .ciTest_mvn_internal(x, statistic=statistic, ...) ## Should <set> go here...
}

## The new one!!
##
.ciTest_mvn_internal <- function(x, statistic="DEV", ...){

  statistic <- match.arg(toupper(statistic), c("DEV","F"))
  if (statistic=="DEV")
    method <- "CHISQ"
  else
    method <- "F"

  S     <- x$cov
  n.obs <- x$n.obs

  vn <- colnames(S)
  K  <- length(vn)

  R   <- vn[-(1:2)]
  v1R <- c(vn[1],R)
  v2R <- c(vn[2],R)

  v1R.idx <- match(v1R, vn)
  v2R.idx <- match(v2R, vn)
  R.idx   <- match(R,   vn)

  d <- n.obs * (log(det(S[v1R.idx, v1R.idx, drop=FALSE])) +
                log(det(S[v2R.idx, v2R.idx, drop=FALSE])) -
                log(det(S[R.idx, R.idx, drop=FALSE])) - log(det(S))
                )

  num.df <- 1
  switch(statistic,
         "DEV"={
           tobs     <- d
           denom.df <- NULL
           p        <- 1-pchisq(tobs, df=num.df)
         },
         "F"={
           tobs     <- (exp(d/n.obs)-1)*(n.obs-K)
           denom.df <- n.obs-K
           p        <- 1-pf(tobs, num.df, denom.df)
         })

  ans <- list(statistic=tobs, p.value=p, df=num.df, denom.df=denom.df,
              statname=statistic, method=method, varNames=vn)
  class(ans) <- "citest"
  ans
}
