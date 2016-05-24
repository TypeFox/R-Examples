internal.1 <- function(semi.major,semi.minor,theta)
{
b.x <- sqrt((semi.major*cos(theta))^2+(semi.minor*sin(theta))^2)
p.a. <- -asin(semi.major*cos(theta)/b.x)+pi/2
retention <- -semi.major*sin(theta)*cos(pi/2+p.a.)+semi.minor*cos(theta)*cos(p.a.)
b.y <- (semi.major*sin(theta)-retention*sin(p.a.))/cos(p.a.)
return(c(b.x,b.y,retention))
}
