#' @rdname f0
#' 
#' @export

iChao1 <- function(f) {
  if(length(f) < 4) f <- c(f, rep(0, 4 - length(f)))
  x <- Chao1(f)
  n <- x["n"]
  s.est <- x["s.obs"] + x["f0"]
  term.1 <- (n - 3) / (4 * n)
  if(f[4] == 0) f[4] <- 1
  term.1 <- term.1 * f[3] / f[4]
  term.2 <- (n - 3) / (2 * (n - 1))
  term.3 <- f[2] * f[3] / f[4]
  term.4 <- f[1] - (term.2 * term.3)
  x["s.est"] <- unname(s.est + term.1 * max(term.4, 0))
  x["f0"] <- unname(x["s.est"] - x["s.obs"])
  x[is.nan(x)] <- NA
  x
}
