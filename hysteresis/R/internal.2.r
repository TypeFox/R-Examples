internal.2 <- function (b.x, b.y, retention,phase.angle) {
p.a. <- atan2(-sqrt((retention^2-b.x^2-b.y^2)^2+4*retention^2*b.y^2)-retention^2+b.x^2+b.y^2,-2*retention*b.y)
  semi.major <- sqrt((b.x*sin(pi/2-p.a.))^2+(b.y*cos(p.a.)+retention*sin(p.a.))^2)
theta <- atan((b.y*cos(p.a.)+retention*sin(p.a.))/(b.x*cos(p.a.)))
theta <- theta + pi*as.numeric(theta < 0)
rotated.angle <- theta*180/pi
semi.minor <- sqrt(b.x^2-b.x^2*sin(pi/2-p.a.)^2)/sin(theta)
return(c(theta,rotated.angle,semi.major,semi.minor))}
