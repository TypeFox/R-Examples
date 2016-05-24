mel1 <- function(cx=32,cy=39,rote.deg=2,semi.major=7,semi.minor=0.23,phase.angle=0,n.points=24,period=24,sd.x=0,sd.y=0) {
theta <- rote.deg/180*pi
 t <-(0:(n.points-1))/period*2*pi
inti <- internal.1(semi.major,semi.minor,abs(theta))
der <- derived.1(semi.major,semi.minor,abs(theta),inti[1],inti[2],inti[3],period)
wrx <- rnorm(n.points,0,sd.x)
wry <- rnorm(n.points,0,sd.y)
x <- semi.major*cos(theta)*cos(t)-semi.minor*sin(theta)*sin(t)+cx+wrx
y <- semi.major*sin(theta)*cos(t)+semi.minor*cos(theta)*sin(t)+cy+wry
ans <- list("values"=c("cx"=cx,"cy"=cy,"rote.deg"=rote.deg,"semi.major"=semi.major,"semi.minor"=semi.minor,"b.x"=inti[1],"b.y"=inti[2],"phase.angle"=phase.angle,
"n.points"=n.points,"period"=period,"area"=der[1],"lag"=der[2],"retention"=inti[3],"coercion"=der[3]),"method"=1,"x"=x,"y"=y)
class(ans) <- "ellipsemake"
ans
}
