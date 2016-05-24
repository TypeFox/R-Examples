nlsboot <-
function(j=NULL,wr1,wr2,x.pred,y.pred,n,cbb,joint){
  if (is.numeric(cbb)==TRUE) {
    xresid2 <- c(wr1,wr1)
    yresid2 <- c(wr2,wr2)
    k <- n/cbb
    xblocks <- sample(1:n,k,replace=TRUE)
    if (joint==FALSE) yblocks <- sample(1:n,k,replace=TRUE)
    else yblocks <- xblocks
    xressamp <- c(t(outer(xblocks,0:(cbb-1),FUN="+")))
    yressamp <- c(t(outer(yblocks,0:(cbb-1),FUN="+")))
    y.boot<-yresid2[yressamp]+y.pred
    x.boot<-xresid2[xressamp]+x.pred
  }
  else {
    if (joint==FALSE) {
    rx <- sample(wr1,n,replace=TRUE)
    ry <- sample(wr2,n,replace=TRUE) 
    }
    else {
      resid.sampler <- sample(1:n,n,replace=TRUE)
      rx <- wr1[resid.sampler]
      ry <- wr2[resid.sampler]
    }
    x.boot<-rx +  x.pred
    y.boot<-ry +  y.pred
  }
          results <- direct(x.boot,y.boot)   
  cx <- as.vector(results$vals["cx"]); cy <- as.vector(results$vals["cy"]); theta <- as.vector(results$vals["theta"]); semi.major <- as.vector(results$vals["semi.major"]); semi.minor <- as.vector(results$vals["semi.minor"]); rote.deg <- as.vector(results$vals["rotated.angle"]);
  
  z <- rep(1,n)   
          nls.polar.fit<-nls(z~
            ( cos(theta)*(x.boot-cx)+sin(theta)*(y.boot-cy))^2/(((ampx)^2+(ampy)^2+((ampx)^2-(ampy)^2)/((cos(theta))^2-(sin(theta))^2))/2)      + (-sin(theta)*(x.boot-cx)+cos(theta)*(y.boot-cy))^2/((ampx^2+ampy^2-(ampx^2-ampy^2)/((cos(theta))^2-(sin(theta))^2))/2),
                             start=list(
                               ampx=sqrt((semi.major*cos(theta))^2+(semi.minor*sin(theta))^2),
                               ampy=sqrt((semi.major*sin(theta))^2+(semi.minor*cos(theta))^2),
                               cx=cx,cy=cy,theta=theta),trace=F
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
          z <- c("cx"=cx,"cy"=cy,"theta"=theta,"semi.major"=semi.major,"semi.minor"=semi.minor,"theta.deg"=rotated.angle)
          z
      }
