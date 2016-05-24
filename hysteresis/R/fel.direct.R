fel.direct <-
function(x,y,ti,period,pred.method){
n <- length(x)
results <- direct(x,y)   
cx <- as.vector(results$vals["cx"]); cy <- as.vector(results$vals["cy"]); theta <- as.vector(results$vals["theta"]); semi.major <- as.vector(results$vals["semi.major"]); semi.minor <- as.vector(results$vals["semi.minor"]); rote.deg <- as.vector(results$vals["rotated.angle"]);
inti <- internal.1(semi.major,semi.minor,theta)
der <- derived.1(semi.major,semi.minor,theta,inti[1],inti[2],inti[3],period)
amps <- derived.amps(inti[1],inti[2],inti[3]) 
focus <- derived.focus(semi.major,semi.minor,theta)

preds <- geometric_distance("x"=x,"y"=y,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=theta,ti,pred.method
                            ,"ampx"=amps[4],"ampy"=amps[5],"lag"=der[2],period)
   
der.summ <- fitstatistics(x,y,preds$pred.x,preds$pred.y,n,method="lm",period)
    
z <- c("cx"=as.vector(cx),"cy"=as.vector(cy),"rote.rad"=as.vector(theta),"semi.major"=as.vector(semi.major),"semi.minor"=as.vector(semi.minor),"rote.deg"=as.vector(rote.deg),
 "area"=der[1],"lag"=der[2],"coercion"=der[3], 
"b.x"=inti[1],"b.y"=inti[2],"retention"=inti[3],"split.angle"=amps[1],"hysteresis.x"=amps[2],"hysteresis.y"=amps[3],"ampx"=amps[4],"ampy"=amps[5], 
        "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3],  "n"=n)

   res <- list("method"="direct","x"=x,"y"=y,"pred.x"=preds$pred.x,"pred.y"=preds$pred.y,"period.time"=preds$period.time,"values"=z,
               "fit.statistics"=der.summ,"residuals"=results$residuals,"fit"=results$fit)
   res
}
