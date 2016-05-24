EMMax.p.HW <-
function(Z,w) {
  p <- apply(Z,2,weighted.mean,w)
  phw <- (4*p[1] + 3*p[2] + 2*p[3] + 1*p[4] + 0*p[5])/4
  c(phw^4,4*phw^3*(1-phw),6*phw^2*(1-phw)^2,4*phw*(1-phw)^3,(1-phw)^4)
}
