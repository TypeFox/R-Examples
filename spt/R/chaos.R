
chaos <- function(abc,...)
  UseMethod("chaos")
chaos.default <- function(abc,...)
  stop("No default function for viewport.  Sorry!")

chaos.spt <- function(abc,iter=10000,...){
  xmin = abc$viewport[1];
  ymin = abc$viewport[2];
  xmax = abc$viewport[3];
  ymax = abc$viewport[4];
  maxA = abc$angles[1]
  main = abc$data.name
  plot(0,0, type='n', bty='n', xaxt='n',yaxt='n',
       xlab='', ylab='', asp=1, main = main,
       xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  polygon(abc$ABC[,1], abc$ABC[,2],lwd=2);
  x0 = (xmax-xmin)*runif(1)+xmin;
  y0 = (ymax-ymin)*runif(1)+ymin;
  if(missing(iter)) iter=10000;
  points(x0,y0,pch='*', col=5)
  points(abc$ABC[,1],abc$ABC[,2],pch=20, col=c(2,3,4))
  ABC = abc$tri
  if(abc$angles[1]>=90) stop("No avilable for obtuse triangle!")
  else .GameSPT(x0,y0,ABC,iter)
  invisible(ABC)
}

chaos.st <- function(abc,iter=10000,...){
  xmin = abc$viewport[1];
  ymin = abc$viewport[2];
  xmax = abc$viewport[3];
  ymax = abc$viewport[4];
  maxA = abc$angles[1]
  main = abc$data.name
  plot(0,0, type='n', bty='n', xaxt='n',yaxt='n',
       xlab='', ylab='', asp=1, main = main,
       xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  polygon(abc$ABC[,1], abc$ABC[,2],lwd=2);
  x0 = (xmax-xmin)*runif(1)+xmin;
  y0 = (ymax-ymin)*runif(1)+ymin;
  if(missing(iter)) iter=10000;
  points(x0,y0,pch='*', col=5)
  points(abc$ABC[,1],abc$ABC[,2],pch=20, col=c(2,3,4))
  ABC = abc$tri
  .GameST(x0,y0,ABC,iter)
  invisible(ABC)
}
