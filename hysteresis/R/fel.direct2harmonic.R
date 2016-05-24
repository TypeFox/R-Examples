fel.direct2harmonic <-
function(x,y,period){
n <- length(x)
results <- direct(x,y)   
cx=as.vector(results$vals["cx"]); cy=as.vector(results$vals["cy"]); theta=as.vector(results$vals["theta"]); semi.major=as.vector(results$vals["semi.major"]); semi.minor=as.vector(results$vals["semi.minor"]); rote.deg=as.vector(results$vals["rotated.angle"]);
inti <- internal.1(semi.major,semi.minor,theta)
der <- derived.1(semi.major,semi.minor,theta,inti[1],inti[2],inti[3],period)
amps <- derived.amps(inti[1],inti[2],inti[3]) 

p.a. <- -asin(semi.major*cos(theta)/inti[1])+pi/2

preds <- geometric_distance("x"=x,"y"=y,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=theta,ti="k",pred.method="find.times"
                            ,"ampx"=amps[4],"ampy"=amps[5],"lag"=der[2],period)
   
ti <- preds$period.time + p.a.
ans2 <- harmonicOptPoint(x,y,ti)
SSE1 <- crossprod(x-preds$pred.x)+crossprod(y-preds$pred.y)
SSE <- ans2$SSE
SSEdiff <- (SSE1-SSE)/SSE

OptCount <- 1
while (SSEdiff > 0.0001) {
  ans2 <- harmonicOptPoint(x,y,ans2$ti)
  SSEdiff <- (SSE - ans2$SSE)/ans2$SSE
  SSE <- ans2$SSE
  OptCount <- OptCount + 1
}

cx <- as.vector(ans2$values["cx"])
cy <- as.vector(ans2$values["cy"])
b.x <- as.vector(ans2$values["b.x"])
b.y <- as.vector(ans2$values["b.y"])
retention <- as.vector(ans2$values["retention"])
phase.angle <- 0
  
inti <- internal.2(b.x,b.y,retention,phase.angle)
n <- length(x)
der <- derived.2(b.x,b.y,retention,period)

residual <- sqrt((x-ans2$pred.x)^2+(y-ans2$pred.y)^2)
amps <- derived.amps(b.x,b.y,retention)  
focus <- derived.focus(inti[3],inti[4],inti[1])

z <- c("cx"=cx,"cy"=cy, "b.x"=b.x,"b.y"=b.y,"phase.angle"=phase.angle, "retention"=retention,
       "area"=der[1], "lag"=der[2], "coercion"=der[3],"rote.rad"=inti[1],
       "rote.deg"=inti[2],"semi.major"=inti[3],"semi.minor"=inti[4],"split.angle"=amps[1],"hysteresis.x"=amps[2],"hysteresis.y"=amps[3],"ampx"=amps[4],"ampy"=amps[5],
       "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3])
res<-list("method"="harmonic2","x"=x,"y"=y,"pred.x"=ans2$pred.x,"pred.y"=ans2$pred.y,"period.time"=ti,
          "values"=z,"fit"=list(ans2,"OptCount"=OptCount,"SSEdiff"=SSEdiff),"fit.statistics"=NULL,"residuals"=residual)
res
}
