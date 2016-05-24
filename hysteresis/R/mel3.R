mel3 <- function(cx=32,cy=39,ampx=6.99,ampy=0.335,lag=2.888,phase.angle=0,n.points=24,period=24,sd.x=0,sd.y=0) {
lag.radian <- lag/period*2*pi
 t <-(0:(n.points-1))/period*2*pi
inti2 <- internal.3(ampx,ampy,lag.radian)
inti1 <- internal.2(inti2[1],inti2[2],inti2[3])
der <- derived.2(inti2[1],inti2[2],inti2[3],period)
wrx <- rnorm(n.points,0,sd.x)
wry <- rnorm(n.points,0,sd.y)
x<-ampx*cos(t+phase.angle)+cx+wrx
y<-ampy*cos(t-lag.radian+phase.angle)+cy+wry
ans <- list("values"=c("cx"=cx,"cy"=cy,"theta.deg"=inti1[2],"semi.major"=inti1[3],"semi.minor"=inti1[4],"b.x"=inti2[1],"b.y"=inti2[2],"phase.angle"=phase.angle,
"n.points"=n.points,"period"=period,"area"=der[1],"lag"=der[2],"retention"=inti2[3],"coercion"=der[3]),"method"=3,"x"=x,"y"=y)
class(ans) <- "ellipsemake"
ans
}
