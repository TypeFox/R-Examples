## rfbaselineScale.R
## compiled in this form: 28/Sept/06


rfbaselineScale <- function(r) {
  dfit <- density(r[r < 3.5*mad(r)])
  modus <- pmin(0,mean(dfit$x[which(dfit$y >= max(dfit$y))]))
  rr <- r[r<= modus] - modus
  rr <- c(rr, -rr)
  sd(rr)
} ## rfbaselineScale
