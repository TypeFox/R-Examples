"SnapPath" <-
  function(step,x, y,f, alpha, alpha0, lambda, Elbow, kernel.function,param.kernel, linear.plot,Size=60,movie=FALSE,movie.root="./",...)
{
  if(movie)jpeg(filename=paste(movie.root,step,".jpg",sep=""),quality=90,width=540,height=540,bg="wheat",...)
  stats<-StatPath(y,f,Elbow)
###only for 2 dim inputs
  x <- x[, 1:2]
  n<-length(y)
    ss<-abs(diff(range(x[,2])))
    soffset<-x*0
    soffset[,2]<-ss/50
    x <- x[, 1:2]
if(!linear.plot){
        xr <- apply(x, 2, range)
      xg <- apply(xr, 2, function(x, Size)
              seq(from = x[1], to = x[2], length = Size), Size)
      xG <- data.matrix(expand.grid(as.list(data.frame(xg))))
      Knew <- kernel.function(x, xG, param.kernel=param.kernel,...)
      fhat <- ((alpha * y) %*% Knew + alpha0)/lambda
      fhat <- matrix(fhat, Size, Size)
      }

    plot(x[,  ], type = "n",xlab="X1",ylab="X2")
   title(paste("              Step: ",format(step,digits=3),"   Error:",format(round(stats$error),digits=3), "   Elbow Size:",format(round(stats$selbow),digits=2),"   Sum Eps:",format(round(stats$margin,2),digits=7)),adj=0)
    pointid <- seq(y)
    points(x[y == 1,  ], col = 2,pch="*",cex=2)
    points(x[y == -1,  ], col = 4,pch="*",cex=2)

    if(n<15)text((x-soffset)[y == 1,  ], labels = paste(pointid[y == 1]), col = 2)
    if(n<15)text((x-soffset)[y == -1,  ], labels = paste(pointid[y == -1]), col = 4)
  
    if(n<15&&length(Elbow))text((x-soffset)[Elbow,  ], labels = paste(pointid[Elbow]), col = 3)
  if(linear.plot){
    beta <- (alpha * y) %*% x
    abline( - alpha0/beta[2],  - beta[1]/beta[2], col = 3, lwd = 3)
    abline(lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
    abline( - lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
  }
  else{
      contour(xg[, 1], xg[, 2], fhat, levels = 0, add = TRUE, labels = "",col=3,lwd=3)
      contour(xg[, 1], xg[, 2], fhat, levels = c(-1,1), add = TRUE, labels = c("",""),col=3,lwd=2,lty=3)
    }
  if(movie)dev.off()
    invisible()
    
}
