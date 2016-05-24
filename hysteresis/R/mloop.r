mloop <- function(cx=0,cy=0,retention=0.2,b.x=0.6,b.y=0.8,n=1,m=1,sd.x=0,sd.y=0,phase.angle=0,n.points=24,period=24,extended.classical=FALSE,seed=NULL) {
 if (!is.null(seed)) set.seed(seed)
 if (extended.classical==FALSE) {
  x<-cx+b.x*cos((1:n.points)/period*2*pi+phase.angle/180*pi)+rnorm(n.points,0,sd.x)
  y<-cy+retention*sin((1:n.points)/period*2*pi+phase.angle/180*pi)^m+b.y*cos((1:n.points)/period*2*pi+phase.angle/180*pi)^n+rnorm(n.points,0,sd.y)
 }
  else {
   direc<-sign(cos((1:n.points)/period*2*pi+phase.angle/180*pi))
   x<-cx+b.x*cos((1:n.points)/period*2*pi+phase.angle/180*pi)+rnorm(n.points,0,sd.x)
   y<-cy+retention*sin((1:n.points)/period*2*pi+phase.angle/180*pi)^m+direc*(b.y*abs(cos((1:n.points)/period*2*pi+phase.angle/180*pi))^n)+rnorm(n.points,0,sd.y)
 }
 if (n==1) beta.split.angle<-atan2(b.y,b.x) 
 else if (n >= 2) beta.split.angle <- 0
  else beta.split.angle<-NA
 hysteresis.x <- 1/sqrt(1+(b.y/retention)^(2/m))
 coercion <- hysteresis.x*b.x
 hysteresis.y <- retention/b.y
  area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs(retention*b.x)
  if ((n%%2)!=1 | (m%%2)!=1) warning("Will not be an actual hysteresis loop if m is not odd, check plot.")
  ans <- list("values"=c("m"=m,"n"=n, "b.x"=b.x,"b.y"=b.y,"phase.angle"=phase.angle,"cx"=cx,"cy"=cy,"retention"=retention,
                         "coercion"=coercion,"area"=area, "beta.split.angle"=beta.split.angle,"hysteresis.x"=hysteresis.x, "hysteresis.y"=hysteresis.y),"x"=x,"y"=y)
class(ans) <- "hysteresisloop"
  ans
}
