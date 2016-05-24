fel.lm <-
function(x,y,ti,period,pred.method){
n <- length(x)
x2 <- -x^2
y2 <- y^2
xy <- x*y
int <- rep(1,n)

model <- lm(terms(x2~0+xy+y2+x+y+int,keep.order=TRUE))
a <- as.vector(c(1,coef(model)))

theta <- atan2(a[2],a[1]-a[3])/2
while(theta<0){theta<-pi/2+theta}


cx <- -(2*a[3]*a[4]-a[2]*a[5])/(4*a[1]*a[3]-a[2]*a[2])
cy <- -(2*a[1]*a[5]-a[2]*a[4])/(4*a[1]*a[3]-a[2]*a[2])

major <- 1/sqrt((a[1]*cos(theta)*cos(theta) + a[2]*cos(theta)*sin(theta) + a[3]*sin(theta)*sin(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))
minor <- 1/sqrt((a[1]*sin(theta)*sin(theta) - a[2]*cos(theta)*sin(theta) + a[3]*cos(theta)*cos(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))


semi.major<- major 
semi.minor<-minor
if (semi.minor > semi.major){
  semi.minor <- semi.major; semi.major <- minor;theta <- theta +pi/2;
}
rotated.angle <- theta*180/pi

inti <- internal.1(semi.major,semi.minor,theta)
der <- derived.1(semi.major,semi.minor,theta,inti[1],inti[2],inti[3],period)
amps <- derived.amps(inti[1],inti[2],inti[3])   
focus <- derived.focus(semi.major,semi.minor,theta)

preds <- geometric_distance("x"=x,"y"=y,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=theta,ti,
                            "pred.method"=pred.method,"ampx"=amps[4],"ampy"=amps[5],"lag"=der[2],period)

der.summ <- fitstatistics(x,y,preds$pred.x,preds$pred.y,n,method="lm",period)
    
z <- c("cx"=cx,"cy"=cy,"rote.rad"=theta,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.deg"=rotated.angle,
 "area"=der[1],"lag"=der[2],"coercion"=der[3],
"b.x"=inti[1],"b.y"=inti[2],"retention"=inti[3],"split.angle"=amps[1],"hysteresis.x"=amps[2],"hysteresis.y"=amps[3],"ampx"=amps[4],"ampy"=amps[5],
         "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3],"n"=n)

   res <- list("method"="lm","x"=x,"y"=y,"pred.x"=preds$pred.x,"pred.y"=preds$pred.y,"period.time"=preds$period.time,"values"=z,
               "fit.statistics"=der.summ,"residuals"=resid(model),"fit"=model)
   res
}
