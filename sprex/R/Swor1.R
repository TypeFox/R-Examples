#' @rdname f0
#' 
#' @export

Swor1 <- function(f, N) {
  if(length(f) == 1) f <- c(f, 0)
  x <- f.stats(f)
  
  # estimate f0
  n <- x["n"]
  q <- n / N
  r <- q / (1 - q)
  w <- n / (n - 1)
  term.1 <- f[1] ^ 2
  term.2 <- 2 * w * f[2]
  term.3 <- r * f[1]
  f0 <- term.1 / (term.2 + term.3)
  if(is.nan(f0)) f0 <- 0

  # est sd(Swor1)
  term.4 <- term.2 * f0 ^ 2
  term.5 <- term.1 * f0
  term.6 <- term.5 ^ 2 / (f[1] ^ 5)
  term.7 <- 4 * w ^ 2 * f[2]
  term.8 <- (f0 / f[1]) ^ 4
  sd.s.est <- unname(sqrt(f0 + term.6 + (term.7 * term.8)))
  if(is.nan(sd.s.est)) sd.s.est <- NA
  
  c(s.est = unname(x["s.obs"] + f0), 
    f0 = unname(f0), 
    x, sd.s.est = unname(sd.s.est)
  )
}
