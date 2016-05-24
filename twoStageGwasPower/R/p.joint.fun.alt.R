p.joint.fun.alt <-
function(xx, c1, c.joint, pi.samples, p0, p1, n.cases, n.controls) {
  # distribution of z.joint, given c1 and c.joint
  fact1 <- termsInt.alt(aa=xx, p0=p0, p1=p1, c.joint=c.joint, pi.samples=pi.samples, n.cases=n.cases, n.controls=n.controls)
  mu1 <- (p1 - p0)/sqrt((p1*(1-p1)/n.cases + p0*(1-p0)/n.controls)/(2*pi.samples))
  
  fact2 <- condnorm(xx, mu1, sd=1, c1=c1)
  result <- fact1$result*fact2
  result
  }
