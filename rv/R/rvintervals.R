

rvintervals <- function (x, rvpoint=rvpar("rvpoint"), ...) {
  which.quantiles <- list(
    "NA" = NA,
    mean = NA,
    median = 0.50,
    "50%" = c(0.25, 0.75),
    "80%" = c(0.10, 0.90),
    "95%" = c(0.025, 0.975)
  )
  .whichq <- function (iv) {
    if (is.numeric(iv)) {
      iv <- paste(100*iv, "%", sep="")
    } else {
      (!is.null(q <- which.quantiles[[iv]])) && return(q)
    }
    if (is.na(iv)) return(NA)
    n <- nchar(iv)
    if (substr(iv,n,n)=="%") {
      ivn <- as.numeric(substr(iv,1,n-1))
      c(0.5-ivn/200, 0.5+ivn/200)
    } else {
      NA
    }
  }
  .length <- function (iv) {
    lg <- .whichq(iv)
    if (length(lg)<=1) 0 else diff(lg)
  }
  probs <- as.vector(na.omit(unlist(lapply(rvpoint, .whichq))))
  if (length(probs)<=1) {
    # A trick to force probs into a named array
    # (won't otherwise return names if we have only one quantile, e.g. 0.50)
    probs <- c(probs, NA)
  }
  if (!all(is.na(probs))) {
    Q <- t(rvquantile(x, probs=probs, ...))
    Q.names <- paste(formatC(100 * probs, format="fg", width=1, digits=3), "%", sep="")
    rownames(Q) <- Q.names
  } else {
    Q <- NULL # DEBUG: will this be ignored if we have only "mean" e.g.? #
  }
  compute.what <- list(
    "NA"   = function () NA,
    mean   = function () t(as.vector(rvmean(x))),
    median = function () Q["50%",,drop=FALSE],
    "50%"  = function () Q[c("25%","75%"),,drop=FALSE],
    "80%"  = function () Q[c("10%","90%"),,drop=FALSE],
    "95%"  = function () Q[c("2.5%","97.5%"),,drop=FALSE]
  )
  .lbl <- function (p) { # From 'quantile.default'
    if (is.null(p) || is.na(p)) return(NA)
    dig <- max(2, getOption("digits"))
    paste(formatC(100 * p, format = "fg", width = 1, digits = dig), "%", sep = "")
  }
  .summaries <- function (iv) {
    if (is.null(f <- compute.what[[iv]])) {
      a <- na.omit(sapply(.whichq(iv),.lbl))
      if (all(a %in% dimnames(Q)[[1]])) {
        return(Q[a,,drop=FALSE])
      } else {
        warning("Cannot understand interval '", iv, "'")
        return(NA)
      }
    }
    a <- f()
    if (is.null(dim(a))) {
       if (length(x)==1) {
          a <- t(a)
       } else {
         na <- names(a)
         dim(a) <- c(1,length(a))
         dimnames(a) <- list(iv, na)
       }
    }
    return(a)
  } 
  lgth <- rev(order(sapply(rvpoint, .length)))
  s <- lapply(rvpoint, .summaries)
  names(s) <- rvpoint
  s[lgth]
}
