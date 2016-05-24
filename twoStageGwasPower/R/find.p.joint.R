find.p.joint <-
function(c.joint, c1, pi.samples, pi.markers) {
  p.joint <- integrate(p.joint.fun, lower=-Inf, upper=-c1, subdivisions=1000, c1=c1, c.joint=c.joint, pi.samples=pi.samples)$value +
           integrate(p.joint.fun, lower=c1, upper=Inf, subdivisions=1000, c1=c1, c.joint=c.joint, pi.samples=pi.samples)$value 
  result <- p.joint
  result
  }
