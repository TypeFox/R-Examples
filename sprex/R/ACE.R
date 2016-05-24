#' @rdname f0
#' 
#' @export

ACE <- function(f) {
  x <- f.stats(f)
  if(length(f) < 10) f <- c(f, rep(0, 10 - length(f)))

  s.rare <- sum(f[1:10])
  x.rare <- sum(1:10 * f[1:10])
  c.ace <- 1 - f[1] / x.rare
  term.1 <- s.rare / c.ace
  term.2 <- sum(1:10 * 0:9 * f[1:10])
  term.3 <- x.rare * (x.rare - 1)
  term.4 <- (term.1 * term.2 / term.3) - 1
  cv.ace <- max(c(term.4, 0))
  f0 <- term.1 + (f[1] * cv.ace / c.ace) - s.rare
  f0 <- f0 * (x["n"] - 1) / x["n"]

  c(s.est = unname(f0 + x["s.obs"]), f0 = unname(f0), x)
}