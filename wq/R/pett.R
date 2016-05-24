pett <- 
function(x, plot = FALSE, order = FALSE, pval = .05, 
         pchs = c(19, 21), ...) {

  # validate args
  if (!is.numeric(x) && !is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a vector, matrix, or data.frame")
  }
  
  # function for single vector
  pet <- function(x) {
  
    # missing data check
    trimna <- cumsum(!is.na(x)) > 0 & rev(cumsum(rev(!is.na(x)))) > 0
    if (is.ts(x)) {
      tx <- time(x)[trimna]
      x1 <- window(x, start = tx[1], end = tx[length(tx)])
    } else {
      x1 <- x[trimna]
    }
    if (anyNA(x1))
      return(setNames(rep(NA, 4), c("pettitt.k", "p.value", 
        "change.point", "change.time")))
    
    # Pettitt change-point statistic
    n <- length(x1)
    d <- sign(outer(x1, x1, "-"))
    u <- sapply(1:(n - 1), function(i) sum(d[1:i, (i + 1):n]))
    pettitt.K <- max(abs(u))
    
    # approximate probability value for Pettitt statistic
    p.value <- 2 * exp(-6 * pettitt.K ^ 2 / (n ^ 3 + n ^ 2))
    p.value <- signif(p.value, 3)
    
    # change position
    change.point <- which.max(abs(u))
    if (is.ts(x1)) {
      change.time <- time(x1)[change.point]
    } else {
      change.time <- change.point
    }
    
    c(
      pettitt.K = pettitt.K, 
      p.value = p.value, 
      change.point = change.point,
      change.time = change.time
    )
  }
  
  # apply pet for each vector
  if (is.null(dim(x))) return(as.list(pet(x)))
  if (identical(ncol(x), 1)) return(as.list(pet(x[, 1])))
  ans <- t(sapply(1:ncol(x), function(i) pet(x[, i])))
  rownames(ans) <- colnames(x)
  
  # plot if TRUE
  if (!plot) {
    ans
  } else {
    if (order) ans <- ans[order(ans[, "change.time"]), ]
    pch <- ifelse(ans[, "p.value"] < pval, pchs[1], pchs[2])
    dotchart(ans[, "change.time"], pch = pch, ...)
  }
}
