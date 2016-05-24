errorplot <-
function(main) {
  plot(0,main=main,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",pch=NA)
  text("not fitted",x=0.5,y=0.5)
}
