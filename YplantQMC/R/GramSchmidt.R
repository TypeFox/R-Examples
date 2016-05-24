#'@author Duncan Murdoch
#'@note Copied from the rgl package, with permission. Will be removed once GramSchmidt exported from rgl in future version.
GramSchmidt <- function (v1, v2, v3, order = 1:3) 
{
  
  xprod <- function (v, w) {
    c(v[2] * w[3] - v[3] * w[2], v[3] * w[1] - v[1] * w[3], v[1] * 
        w[2] - v[2] * w[1])
  }
  veclen <- function (v)sqrt(sum(v^2))
  normalize <- function (v)v/veclen(v)
  
  A <- rbind(v1, v2, v3)
  A <- A[order, ]
  v1 <- A[1, ]
  v2 <- A[2, ]
  v3 <- A[3, ]
  if (isTRUE(all.equal(as.numeric(v1), c(0, 0, 0)))) 
    v1 <- xprod(v2, v3)
  v1 <- normalize(v1)
  v2 <- v2 - sum(v2 * v1) * v1
  if (isTRUE(all.equal(as.numeric(v2), c(0, 0, 0)))) 
    v2 <- xprod(v3, v1)
  v2 <- normalize(v2)
  v3 <- v3 - sum(v3 * v1) * v1 - sum(v3 * v2) * v2
  if (isTRUE(all.equal(as.numeric(v3), c(0, 0, 0)))) 
    v3 <- xprod(v1, v2)
  v3 <- normalize(v3)
  rbind(v1, v2, v3)[order(order), ]
}
