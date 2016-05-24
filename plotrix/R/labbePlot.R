# labbePlot reads a list of comparison trial results between two conditions,
# typically "placebo" and "intervention" for clinical trials,
# and plots each result as a circle centered at the percentage of successes
# for each condition, with radius proportional to the square root of the
# number of observations (multiplied by "circle.mag" for very small or large
# numbers of observations), with optional fill colors.

labbePlot<-function(x,main="L'Abbe plot",
 xlab="Percent positive response with placebo",
 ylab="Percent positive response with treatment",
 labels=NULL,col=NA,circle.mag=0.5,add=FALSE,...) {

 if(is.list(x)) {
  if(!add)
   plot(0,xlim=c(0,100),ylim=c(0,100),main=main,xlab=xlab,
    ylab=ylab,type="n",...)
  for(trial in 1:length(x)) {
   if(is.matrix(x[[trial]])) {
    sum_treat<-sum(x[[trial]][1,])
    sum_interv<-sum(x[[trial]][2,])
    xpos<-100*x[[trial]][1,1]/sum_treat
    ypos<-100*x[[trial]][2,1]/sum_interv
    rad<-circle.mag*sqrt(sum_treat+sum_interv)
   }
   else {
    xpos<-x[[trial]][1]
    ypos<-x[[trial]][2]
    rad<-circle.mag*sqrt(x[[trial]][3])
   }
   circle.col<-ifelse(is.list(col),col[[trial]],col)
   draw.circle(xpos,ypos,rad,col=circle.col)
   if(!is.null(labels[[trial]])) {
    textcol<-ifelse(colSums(col2rgb(circle.col)*c(1,1.4,0.6)) < 350,
     "white", "black")
    text(xpos,ypos,labels[[trial]],col=textcol)
   }
  }
  segments(0,0,100,100)
 }
 else {
  cat("labbePlot: x must be a list of 2x2 tables OR\n")
  cat("3 element numeric vectors of percent, percent, N (see help page)\n")
 }
}
