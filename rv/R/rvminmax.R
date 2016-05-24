
rvminmax <- function (x, na.rm=FALSE, finite=FALSE, defaults=c(-Inf,Inf)) {
  ## NOEXPORT
  .minmax <- function (x, FUN, def) {
    if (na.rm) x <- x[!is.na(x)]
    if (finite) x <- x[is.finite(x)]
    if (length(x)==0) def else FUN(x)
  }
  m0 <- rvsimapply(x, .minmax, min, defaults[1])
  m1 <- rvsimapply(x, .minmax, max, defaults[2])
  list(min=m0, max=m1)
}

