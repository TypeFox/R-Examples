# x is decreasing
# x would be lambda
interp_w <- function(s, x) {
  n <- length(x)
  s_le_x <- s <= x # s <= x?
  if (!any(s_le_x)) return(list("indices"=c(1L, 1L), "weights"=c(0, 1)))
  lower <- max(which(s_le_x))
  if (lower == n) return(list("indices"=c(n, n), "weights"=c(0, 1)))
  # so  1 <= lower < n
  upper <- lower + 1L
  upper_w <- x[lower]-s
  lower_w <- s-x[upper]
  return(list("indices"=c(lower, upper),
              "weights"=c(lower_w, upper_w) / (lower_w + upper_w)))
}


"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y