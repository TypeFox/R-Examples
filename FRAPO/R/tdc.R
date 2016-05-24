##
## Tail dependence coefficient
## (non-parametric estimation)
##
tdc <- function(x, method = c("EmpTC", "EVT"), lower = TRUE, k = NULL, ...){
  x <- as.matrix(x)
  ifelse(is.null(k), k <- floor(sqrt(nrow(x))), k <- as.integer(k))
  method = match.arg(method)
  m <- nrow(x)
  N <- ncol(x)
  idx <- combn(1:N, 2)
  r <- apply(x, 2, rank, ...)
  if(method == "EmpTC"){
    if(lower){
      td <- apply(idx, 2, function(y)
                  sum((r[, y[1]] <= k) & (r[, y[2]] <= k)) / k)
    } else {
      td <- apply(idx, 2, function(y)
                  sum((r[, y[1]] > m - k) & (r[, y[2]] > m - k)) / k)
    }
  }
  if(method == "EVT"){
    if(lower){
      td <- apply(idx, 2, function(y)
                  2 - sum((r[, y[1]] <= k) | (r[, y[2]] <= k)) / k)
    } else {
      td <- apply(idx, 2, function(y)
                  2 - sum((r[, y[1]] > m - k) | (r[, y[2]] > m - k)) / k)
    }
  } 
  tdm <- diag(N)
  tdm[t(idx)] <- td
  tdm[lower.tri(tdm)] <- td
  colnames(tdm) <- rownames(tdm) <- colnames(x)
  return(tdm)
}
