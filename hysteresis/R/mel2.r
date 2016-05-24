mel2 <- function(cx=32,cy=39,b.x=6.99,b.y=0.244,retention=0.23,phase.angle=0,n.points=24,period=24,sd.x=0,sd.y=0) {
ti <- (0:(n.points-1))*2/period*pi
wrx <- rnorm(n.points,0,sd.x)
wry <- rnorm(n.points,0,sd.y)
x <- b.x*cos(ti+phase.angle)+cx+wrx
y <- b.y*cos(ti+phase.angle)+retention*sin(ti+phase.angle)+cy+wry
inti <- internal.2(b.x,b.y,retention,phase.angle)
der <- derived.2(b.x,b.y,retention,period)
ans <- list("values"=c("cx"=cx,"cy"=cy,"rote.deg"=inti[2],"semi.major"=inti[3],"semi.minor"=inti[4],"b.x"=b.x,"b.y"=b.y,"phase.angle"=phase.angle,
"n.points"=n.points,"period"=period,"area"=der[1],"lag"=der[2],"retention"=retention,"coercion"=der[3]),"method"=2,"x"=x,"y"=y)
class(ans) <- "ellipsemake"
ans
}
