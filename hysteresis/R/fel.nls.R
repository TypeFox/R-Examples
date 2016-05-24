fel.nls <-
function (x,y,ti,control,period,pred.method) {
n <- length(x)
results <- direct(x,y)   
cx <- as.vector(results$vals["cx"]); cy <- as.vector(results$vals["cy"]); theta <- as.vector(results$vals["theta"]); semi.major <- as.vector(results$vals["semi.major"]); semi.minor <- as.vector(results$vals["semi.minor"]); rote.deg <- as.vector(results$vals["rotated.angle"]);
z <- rep(1,n)   
nls.polar.fit<-nls(z~
     ( cos(theta)*(x-cx)+sin(theta)*(y-cy))^2/(((ampx)^2+(ampy)^2+((ampx)^2-(ampy)^2)/((cos(theta))^2-(sin(theta))^2))/2)      + (-sin(theta)*(x-cx)+cos(theta)*(y-cy))^2/((ampx^2+ampy^2-(ampx^2-ampy^2)/((cos(theta))^2-(sin(theta))^2))/2),
     start=list(
     ampx=sqrt((semi.major*cos(theta))^2+(semi.minor*sin(theta))^2),
     ampy=sqrt((semi.major*sin(theta))^2+(semi.minor*cos(theta))^2),
     cx=cx,cy=cy,theta=theta),control=control,trace=F
    )

    #summary(nls.polar.fit)
    ampx<-as.vector(coef(nls.polar.fit)[1])
    ampy<-as.vector(coef(nls.polar.fit)[2])

    theta<-as.vector(coef(nls.polar.fit)[5])
    rotated.angle <- theta*180/pi
    
    cx <- as.vector(coef(nls.polar.fit)[3])
    cy <- as.vector(coef(nls.polar.fit)[4])
                                                                  
    semi.major<-sqrt((((ampx)^2+(ampy)^2+((ampx)^2-(ampy)^2)/((cos(theta))^2-(sin(theta))^2))/2)  )
    semi.minor<-sqrt((ampx^2+ampy^2-(ampx^2-ampy^2)/((cos(theta))^2-(sin(theta))^2))/2)


if (semi.minor > semi.major){
  lam <- semi.minor; semi.minor <- semi.major; semi.major <- lam;theta <- theta +pi/2;rotated.angle <- theta/180*pi;
}
    MSE<-(summary(nls.polar.fit)$sigma)^2
    
    
inti <- internal.1(semi.major,semi.minor,theta)
der <- derived.1(semi.major,semi.minor,theta,inti[1],inti[2],inti[3],period)
split.angle <- atan2(sqrt(ampy^2-inti[3]^2),ampx)
split.angle <- split.angle*180/pi
hysteresis.y <- inti[3]/inti[2]
hysteresis.x <- 1/sqrt(1+(inti[2]/inti[3])^2)
focus <- derived.focus(semi.major,semi.minor,theta)
preds <- geometric_distance("x"=x,"y"=y,"cx"=cx,"cy"=cy,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.rad"=theta,ti,"pred.method"=pred.method,"ampx"=ampx,"ampy"=ampy,"lag"=der[2],"period"=period)    
der.summ <- fitstatistics(x,y,preds$pred.x,preds$pred.y,n,method="nls",period)
    
z <- c("cx"=cx,"cy"=cy,"rote.rad"=theta,"semi.major"=semi.major,"semi.minor"=semi.minor,"rote.deg"=rotated.angle,
"area"=der[1],"lag"=der[2],"coercion"=der[3], 
        "b.x"=inti[1],"b.y"=inti[2],"retention"=inti[3],"split.angle"=split.angle,"hysteresis.x"=hysteresis.x,"hysteresis.y"=hysteresis.y, "ampx"=ampx,"ampy"=ampy, 
        "focus.x"=focus[1],"focus.y"=focus[2], "eccentricity"=focus[3],"n"=n)
   res<-list("method"="nls","x"=x,"y"=y,"pred.x"=preds$pred.x,"pred.y"=preds$pred.y,"period.time"=preds$period.time,
             "values"=z,"fit.statistics"=c("MSE.nls"=MSE,der.summ),
             "residuals"=resid(nls.polar.fit),"fit"=nls.polar.fit)
   res
   }
