lengthKey<-function(x,y,tickpos,scale) {
 par(xpd=TRUE)
 nticks<-length(tickpos)
 segments(x,y,x+tickpos[nticks]*scale,y)
 for(tick in 1:nticks) {
  segments(x+tickpos[tick]*scale,y,
   x+tickpos[tick]*scale,y+tickpos[nticks]*scale/20)
  text(x+tickpos[tick]*scale,y+tickpos[nticks]*scale/15,
  tickpos[tick],adj=c(0.5,0))
 }
 par(xpd=FALSE)
}
