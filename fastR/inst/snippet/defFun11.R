x = (0:10)/10

myData=data.frame(x=x, y=sin(x))

panel.xyplotWithDiag <- function(x,y,...) {
    panel.xyplot(x,y,...)
    panel.abline(a=0,b=1,col="gray30",lwd=2)
}
