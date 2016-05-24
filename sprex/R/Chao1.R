#' @rdname f0
#' 
#' @export

Chao1 <- function(f) {
  x <- f.stats(f)
  s.obs <- unname(x["s.obs"])
  if(length(f) == 1) f <- c(f, 0)
  
  f0 <- if(f[2] > 0) {
    f[1] ^ 2 / (2 * f[2])
  } else {
    term.1 <- f[1] * (f[1] - 1)
    term.2 <- 2 * (f[2] + 1)
    term.1 / term.2
  }
  f0 <- f0 * (x["n"] - 1) / x["n"]
  
  x <- c(s.est = unname(f0 + s.obs), f0 = unname(f0), x)
  x[is.nan(x)] <- NA
  x
}