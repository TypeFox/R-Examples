"precip.dep" <-
function(x, norm, var = "precip") {
  orig <- as.character(substitute(x))[[1]]
  sc <- seas.df.check(x, orig, var)
  if (!inherits(norm, "seas.norm"))
    stop(gettextf("%s is not a %s object",
                  sQuote(orig), sQuote("seas.norm")))
  x$fact <- mkseas(x, norm$width)
  x$dep <- x[[var]]
  na <- is.na(x[[var]])
  x$dep[na] <- 0
  for (i in levels(x$fact)) {
    sl <- x$fact == i & !na
    n <- norm$seas[i, var]
    x$dep[sl] <- x$dep[sl] - n
  }
  x$dep <- cumsum(x$dep)
  x
}
