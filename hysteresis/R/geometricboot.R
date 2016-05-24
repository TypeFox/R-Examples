geometricboot <-
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
  x <- x.boot
  y <- y.boot
  start <- direct(x,y) 
  
  ti<-n
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
  results <- geom_ellipse(model,1.001)  
  cx <- as.vector(results$values[4]); cy <- as.vector(results$values[5]); 
  theta <- as.vector(results$values[1]); semi.major <- as.vector(results$values[2]); 
  semi.minor <- as.vector(results$values[3]); 
          z <- c("cx"=cx,"cy"=cy,"theta"=theta,"semi.major"=semi.major,"semi.minor"=semi.minor,"theta.deg"=theta*180/pi,"phase.angle"=ti[1])
          z
      }
