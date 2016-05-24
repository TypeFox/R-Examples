fel.geometric <-
function(x,y,control,period){
start <- direct(x,y) 

ti<-numeric(length(x))
for (i in 1:length(x)) {
  x0<-x[i]
  y0<-y[i]
  zmin1<-optimize(ellipsespot,c(0,pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  zmin2<-optimize(ellipsespot,c(pi,2*pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  ti[i]<-ifelse(zmin1$objective < zmin2$objective, zmin1, zmin2)[[1]]
}
pred.x<-start$vals["cx"] +start$vals["semi.major"]*cos(start$vals["theta"])*cos(ti)-start$vals["semi.minor"]*sin(start$vals["theta"])*sin(ti)
pred.y<-start$vals["cy"] +start$vals["semi.major"]*sin(start$vals["theta"])*cos(ti)+start$vals["semi.minor"]*cos(start$vals["theta"])*sin(ti)

model <- list("period.time"=ti,"values"=c("cx"=as.vector(start$vals["cx"]),"cy"=as.vector(start$vals["cy"]),
                                    "semi.major"=as.vector(start$vals["semi.major"]),"semi.minor"=as.vector(start$vals["semi.minor"]),
                                             "rote.rad"=as.vector(start$vals["theta"])),"x"=x,"y"=y)
results <- geom_ellipse(model,control)  
n <- length(model$x)
cx <- as.vector(results$values[4]); cy <- as.vector(results$values[5]); 
theta <- as.vector(results$values[1]); semi.major <- as.vector(results$values[2]); 
semi.minor <- as.vector(results$values[3]); 
if (theta < 0) while (theta < 0) theta <- theta + 2*pi
if (theta > 2*pi) while (theta > 2*pi) theta <- theta - 2*pi
inti <- internal.1(semi.major,semi.minor,theta)
der <- derived.1(semi.major,semi.minor,theta,inti[1],inti[2],inti[3],period)
amps <- derived.amps(inti[1],inti[2],inti[3]) 
focus <- derived.focus(semi.major,semi.minor,theta)

pred.x <- model$x-results[[3]]
pred.y <- model$y - results[[4]]
   
der.summ <- fitstatistics(model$x,model$y,pred.x,pred.y,n,method="lm",period)
    
z <- c("cx"=as.vector(cx),"cy"=as.vector(cy),"rote.rad"=as.vector(theta),"semi.major"=as.vector(semi.major),"semi.minor"=semi.minor,"rote.deg"=as.vector(theta)*180/pi,
 "phase.angle"=ti[1] ,"area"=der[1],"lag"=der[2],"coercion"=der[3], 
"b.x"=inti[1],"b.y"=inti[2],"retention"=inti[3],"split.angle"=amps[1],"hysteresis.x"=amps[2],"hysteresis.y"=amps[3],"ampx"=amps[4],"ampy"=amps[5], 
        "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3],  "n"=n)

   res <- list("method"="geometric","x"=x,"y"=y,"pred.x"=pred.x,"pred.y"=pred.y,"period.time"=results$period.time,"values"=z,
               "fit.statistics"=der.summ,"residuals"=sqrt(results[[3]]^2+results[[4]]^2),"fit"=results)

class(res) <- "ellipsefit"  
res
}
