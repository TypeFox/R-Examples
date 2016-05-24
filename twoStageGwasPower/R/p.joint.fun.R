p.joint.fun <-
function(xx, c1, c.joint, pi.samples) {
  # distribution of z.joint, given c1 and c.joint
  fact1 <- termsInt(aa=xx, c.joint=c.joint, pi.samples)
  fact2 <- condnorm(xx, c1=c1)
  result <- fact1$result*fact2
  result
  }
