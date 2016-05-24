mloop2r <- function(cx=0,cy=0,retention.above=0.2,retention.below=0.15,b.x=0.6,b.y=0.8,n=1,m=1,sd.x=0,sd.y=0,phase.angle=0,n.points=24,period=24,extended.classical=FALSE,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  t2 <- (1:n.points)/period*2*pi+phase.angle/180*pi
  Ind <- (t2 < pi) & (t2 > 0) 
  if (extended.classical==FALSE) {
  x<-cx+b.x*cos(t2)+rnorm(n.points,0,sd.x)
  y<-cy+retention.above*Ind*sin(t2)^m+retention.below*(1-Ind)*sin(t2)^m+b.y*cos((1:n.points)/period*2*pi+phase.angle/180*pi)^n+rnorm(n.points,0,sd.y)
 }
  else {
   direc<-sign(cos(t2))
   x<-cx+b.x*cos(t2)+rnorm(n.points,0,sd.x)
   y<-cy+retention.above*Ind*sin(t2)^m+retention.below*(1-Ind)*sin(t2)^m+direc*(b.y*abs(cos(t2))^n)+rnorm(n.points,0,sd.y)
 }
  if (n==1) beta.split.angle<-atan2(b.y,b.x)*180/pi 
  else if (n >= 2) beta.split.angle <- 0
  else beta.split.angle<-NA
  hysteresis.x.above <- 1/sqrt(1+(b.y/retention.above)^(2/m))
  hysteresis.x.below <- 1/sqrt(1+(b.y/retention.below)^(2/m))
  coercion.above <- hysteresis.x.above*b.x
  coercion.below <- hysteresis.x.below*b.x
  hysteresis.y.above <- retention.above/b.y
  hysteresis.y.below <- retention.below/b.y
  area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs((retention.above+retention.below)*b.x)/2
  lag.above<-abs(atan2(retention.above,b.y))*period/(pi*2)
  lag.below<-abs(atan2(retention.below,b.y))*period/(pi*2)
  if ((n%%2)!=1 | (m%%2)!=1) warning("Will not be an actual hysteresis loop if m is not odd, check plot.")
  ans <- list("values"=c("n"=n, "m"=m,"b.x"=b.x,"b.y"=b.y,"phase.angle"=as.vector(phase.angle),"cx"=cx,"cy"=cy,"retention.above"=retention.above,
                         "retention.below"=retention.below, "coercion.above"=coercion.above,"coercion.below"=coercion.below,"area"=area, "lag.above"=lag.above,"lag.below"=lag.below,"beta.split.angle"=beta.split.angle,
                         "hysteresis.x.above"=hysteresis.x.above,"hysteresis.x.below"=hysteresis.x.below, "hysteresis.y.above"=hysteresis.y.above,"hysteresis.y.below"=hysteresis.y.below),"x"=x,"y"=y)
class(ans) <- "hysteresisloop"
  ans
}
