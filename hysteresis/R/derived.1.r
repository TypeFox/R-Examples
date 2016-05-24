derived.1 <- function(semi.major,semi.minor,theta,b.x,b.y,retention,period) {
area <- semi.major*semi.minor*pi
coercion <- b.x/(sqrt(1+(b.y/retention)^2))
lag<-abs(atan2(retention,b.y))*period/(pi*2)
z <- c(area,lag,coercion)
z
}
