find.p.joint.alt <-
function(c.joint, c1, pi.samples, pi.markers, p0, p1, n.cases, n.controls) {
  
  p.joint <- integrate(p.joint.fun.alt, lower=-Inf, upper=-c1, subdivisions=1000, c1=c1, c.joint=c.joint, 
                    pi.samples=pi.samples, p0=p0, p1=p1, n.cases=n.cases, n.controls=n.controls)$value +
           integrate(p.joint.fun.alt, lower=c1, upper=Inf, subdivisions=1000, c1=c1, c.joint=c.joint, 
                     pi.samples=pi.samples, p0=p0, p1=p1, n.cases=n.cases, n.controls=n.controls)$value 
  result <- p.joint
  result
  }
