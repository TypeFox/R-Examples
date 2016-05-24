mel4 <- function(x2=0.002293,xy=-.06960,y2=0.9976,x=2.567,y=-75.58,int=1432.7,
                 phase.angle=0,n.points=24,period=24,sd.x=0,sd.y=0) {
  a<-c(x2,xy,y2,x,y,int)
  theta <- atan2(a[2],a[1]-a[3])/2
  while(theta<0){theta<-pi/2+theta} 
  
  a[4] <- abs(a[4])
  a[5] <- abs(a[5])
  
  cx <- -(2*a[3]*a[4]-a[2]*a[5])/(4*a[1]*a[3]-a[2]*a[2])
  cy <- -(2*a[1]*a[5]-a[2]*a[4])/(4*a[1]*a[3]-a[2]*a[2])
  
  major <- 1/sqrt((a[1]*cos(theta)*cos(theta) + a[2]*cos(theta)*sin(theta) + a[3]*sin(theta)*sin(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))
  minor <- 1/sqrt((a[1]*sin(theta)*sin(theta) - a[2]*cos(theta)*sin(theta) + a[3]*cos(theta)*cos(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))
  
  
  semi.major<- abs(major) 
  semi.minor<-abs(minor)
  if (semi.minor > semi.major){
    semi.minor <- semi.major; semi.major <- abs(minor);theta <- theta +pi/2;
  }
   rotated.angle <-theta*180/pi
  t <-(0:(n.points-1))/period*2*pi
  inti <- internal.1(semi.major,semi.minor,abs(theta))
  der <- derived.1(semi.major,semi.minor,abs(theta),inti[1],inti[2],inti[3],period)
  wrx <- rnorm(n.points,0,sd.x)
  wry <- rnorm(n.points,0,sd.y)
  x <- semi.major*cos(theta)*cos(t)-semi.minor*sin(theta)*sin(t)+cx+wrx
  y <- semi.major*sin(theta)*cos(t)+semi.minor*cos(theta)*sin(t)+cy+wry
  ans <- list("values"=c("cx"=cx,"cy"=cy,"rote.deg"=rotated.angle,"semi.major"=semi.major,"semi.minor"=semi.minor,"b.x"=inti[1],"b.y"=inti[2],"phase.angle"=phase.angle,
                         "n.points"=n.points,"period"=period,"area"=der[1],"lag"=der[2],"retention"=inti[3],"coercion"=der[3]),"method"=1,"x"=x,"y"=y)
  class(ans) <- "ellipsemake"
  ans
}
