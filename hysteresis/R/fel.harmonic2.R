fel.harmonic2 <-
function(x,y,ti,period,pred.method){


  Ta.lm<-lm.fit(cbind(rep(1,length(x)),sin(ti),cos(ti)),x)            
  b.x<-sqrt(coef(Ta.lm)[[2]]^2+coef(Ta.lm)[[3]]^2)
  phase.angle<- atan2(coef(Ta.lm)[[3]],coef(Ta.lm)[[2]])-pi/2
  rad<-ti+phase.angle
  cx<-coef(Ta.lm)[[1]]
  Tb.lm<-lm.fit(cbind(rep(1,length(y)),sin(rad),cos(rad)),y)

  b.y<-coef(Tb.lm)[[3]]
  retention<- coef(Tb.lm)[[2]]
  cy<-coef(Tb.lm)[[1]]
inti <- internal.2(b.x,b.y,retention,phase.angle)
  n <- length(x)
 der <- derived.2(b.x,b.y,retention,period)
    
     newx<-b.x*cos(rad)+cx
 
newy<-b.y*cos(rad)+retention*sin(rad)+cy

der.summ <- fitstatistics(x,y,newx,newy,n,method="harmonic2",period)
  residual <- sqrt((x-newx)^2+(y-newy)^2)
amps <- derived.amps(b.x,b.y,retention)  
  focus <- derived.focus(inti[3],inti[4],inti[1])
  
   z <- c("cx"=cx,"cy"=cy, "b.x"=b.x,"b.y"=b.y,"phase.angle"=phase.angle, "retention"=retention,
          "area"=der[1], "lag"=der[2], "coercion"=der[3],"rote.rad"=inti[1],
          "rote.deg"=inti[2],"semi.major"=inti[3],"semi.minor"=inti[4],"split.angle"=amps[1],"hysteresis.x"=amps[2],"hysteresis.y"=amps[3],"ampx"=amps[4],"ampy"=amps[5],
          "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3])
   res<-list("method"="harmonic2","x"=x,"y"=y,"pred.x"=newx,"pred.y"=newy,"period.time"=ti,
             "values"=z,"fit"=list(Ta.lm,Tb.lm),"fit.statistics"=der.summ,"residuals"=residual)
   res
    }
