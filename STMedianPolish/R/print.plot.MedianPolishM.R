#' @export 
#' 
print.plot.MedianPolishM <-
function(x,...)
{

if(x$Graphic != 2)
{
layout(matrix(c(1:length(dim(x$residuals))), 1, length(dim(x$residuals)), byrow = TRUE))
a<-length(dimnames(x$residuals))

if(a!=0)
{
for(i in 1:length(dim(x$residuals))){
stripchart(t(x[[i]])~as.factor(dimnames(x$residuals)[[i]]),
,xlab="Factors", ylab=names(dimnames(x$residuals))[i],pch=i,vertical=T)}
}else
{
for(i in 1:1:length(dim(x$residuals))){
stripchart(t(x[[i]])~as.factor(1:dim(x$residuals)[i]),
,xlab="Factors", ylab=names(dimnames(x$residuals))[i],pch=i,vertical=T)}
}
}else
{
layout(matrix(c(1,2,3,4,4,5), 2, 3, byrow = TRUE))

COLORS <-c(rep("black",x$iter-1),"red")
LTY <-c(rep(2,x$iter-1),1)
Thick <-c(rep(1,x$iter-1),2)

TEF<-c("x","y","z","Time")

for(i in 1:4){
ts.plot(t(x[[i]]),xlab="Intervals", ylab="Effect",
main=paste(TEF[i]," effect"),col=COLORS,lty=LTY,lwd=Thick)
}
plot(1, type = "n", axes = FALSE,xlab="",ylab="")
legend("bottom", cex=1.2,pt.cex=3,title="Legend",  bty = "n",
       legend=c( "Iterations","Convergence"),col=c("black","red"), lty=c(2,1),lwd=c(1,2))
}
}
