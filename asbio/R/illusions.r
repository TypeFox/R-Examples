illusions<-function(ill.no=1) {
old.par <- par(no.readonly = TRUE)
if(ill.no==1){
#Illusion#1
par(mar = c(5,0,4,0))
plot(seq(1,10),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",main="Which line is longer?")
arrows(2.5,4,8.5,4,lwd=4)
arrows(8.5,4,2.5,4,lwd=4)
arrows(8.5,6,2.5,6,angle=140,lwd=4)
arrows(2.5,6,8.5,6,angle=140,lwd=4)
readline("Press return for next plot")
rect(1,1,2.5,10, col = "white", border = "white")
rect(8.5,1,10,10, col = "white", border = "white")
rect(2.5,1,8.5,5, col = "white", border = "white")
segments(2.5,4,8.5,4,lwd=4)
abline(v=c(2.5,8.5), col = "gray", lty = 2)
}
if(ill.no==2){
#Illusion#2
plot(seq(1,10),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",main="Are the red lines straight?")
points(5,5,pch=8,cex=47,lwd=8)
segments(5,5,3,0,lwd=8)
segments(5,5,7,0,lwd=8)
segments(5,5,0,3,lwd=8)
segments(5,5,0,7,lwd=8)
segments(5,5,3,10,lwd=8)
segments(5,5,7,10,lwd=8)
segments(5,5,10,3,lwd=8)
segments(5,5,10,7,lwd=8)
segments(c(2,4,6,8),c(1,1,1,1),c(2,4,6,8),c(10,10,10,10),lwd=8,col="red")
}
if(ill.no==3){
#Illusion#3
op<-par(mar = c(0,0,0,0), bg = "black")
pn<-plot.new()
x = seq(0, 1, length = 10)
y = seq(0, 1, length = 10)
abline(v = x, h = y, col = "gray", lwd = 6)
points(rep(x, each=10), rep(y, 10), col = "white", 
cex = 3, pch = 20)
par(op)
print("Code for this illusion follow from Yihui Xie's package animation")
}
on.exit(par(old.par))
}
