seaKen <- 
function(x, plot = FALSE, type = c("slope", "relative"), order = FALSE, 
         pval = .05, mval = .5, pchs = c(19, 21), ...) {

  # validate args
  if (!is.numeric(x) && !is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a vector, matrix, or data.frame")
  if (!is.null(ncol(x)) && is.null(colnames(x)))
    colnames(x) <- paste("series_", 1:ncol(x), sep="")
  type <- match.arg(type)

  # test for single series
  sk <- function(x) {
    
    # validate args
    if (!is.ts(x))
      stop("'x' must be of class 'ts'")
    if (identical(frequency(x), 1))
      stop("'x' must be a seasonal time series with frequency > 1")
    
    # extend series to full years
    fr <- frequency(x)
    xmod <- length(x) %% fr
    if (!identical(xmod, 0))
      x <- ts(c(x, rep(NA, fr - xmod)), start = start(x), frequency = fr)

    # apply mannKen to matrix of months
    x1 <- matrix(x, ncol = fr, byrow = TRUE)
    mk1 <- mannKen(x1)
  
    # calculate sen slope
    sen.slope <- median(mk1[, "sen.slope"], na.rm = TRUE)
    sen.slope.rel <- sen.slope / abs(median(x, na.rm = TRUE))
    
    # calculate sen slope significance
    S <- sum(mk1[, "S"])
    varS <- sum(mk1[, "varS"])
    Z <- (S - sign(S)) / sqrt(varS)
    p.value <- 2 * pnorm(-abs(Z))
    
    miss <- round(sum(mk1[, "miss"] >= .5) / fr, 3)
    c(sen.slope = sen.slope, 
      sen.slope.rel = sen.slope.rel,
      p.value = p.value, 
      miss = miss)
  }
	
  # apply sk for each series
  if (is.null(dim(x))) return(as.list(sk(x)))
  if (ncol(x) == 1) return(as.list(sk(x[, 1])))
  ans <- t(sapply(1:ncol(x), function(i) sk(x[, i])))
  rownames(ans) <- colnames(x)
  
	# plot if TRUE
	if (!plot) {
	  ans
	} else {
	  v1 <- switch(type,
	               slope = "sen.slope",
	               relative = "sen.slope.rel")
	  if (order) ans <- ans[order(ans[, v1]), ]
	  pch <- ifelse(ans[, "miss"] >= mval, NA, 
	                ifelse(ans[, "p.value"] < pval, pchs[1], pchs[2]))
	  dotchart(ans[, v1], pch = pch, ...)
	}
}
