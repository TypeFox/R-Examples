local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# R code for the tranform example: sine of a normal variate

plotSN <- function() {
  getWinVal(scope="L");
  x1 <- rnorm(n,mu,sig); x2 <- rnorm(n,mu,sig); xx <- c(x1,x2);
  y1 <- sin(2*pi*x1); y2 <- sin(2*pi*x2); yy <- c(y1,y2);
  mx <- mean(xx); sx <- sd(xx);
  my <- mean(yy); sy <- sd(yy);
  xy <- c(mx,my,sx,sy); xy <- round(1000*xy)/1000;
  setWinVal(list(mx=xy[1],my=xy[2],sx=xy[3],sy=xy[4]));
  xlo <- min(0,xx); xhi <- max(1,xx);
  x <- seq(xlo,xhi,length=1000); y <- sin(2*pi*x);
  xpars <- paste("(xm = ",mu,", xs = ",sig,")",sep="");
  resetGraph();
  if( ptype=="c" ) {
    plot(x,y,type="l",xlab=paste("x",xpars));
    dy <- 0.05;
    yu <- runif(2*n,-dy,dy);
    points(xx,yy,col="red");
    points(xx,yu,cex=.5,col="blue"); };
  if( ptype=="p" ) {
    plot(y1,y2,xlab=paste("y1",xpars),pch=19,cex=.1); };
  if( ptype=="h" ) {
    hist(yy,xlab=paste("y",xpars),main=""); };
  };

require(PBSmodelling); createWin("SineNormWin.txt")

}) # end local scope
