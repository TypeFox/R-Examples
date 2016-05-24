find.p.joint.alpha <-
function(c.joint, c1, pi.samples, pi.markers, alpha.marker) {
  result.A <- find.p.joint(c.joint=c.joint, c1=c1, pi.samples=pi.samples, pi.markers=pi.markers)
  #M <- 2.380*10^6
  #result <- (result.A - alpha.marker/(M*pi.markers))^2
  result <- (result.A - alpha.marker/(pi.markers))^2
  result
  }
