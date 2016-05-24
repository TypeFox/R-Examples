xyplot(y~x, data=myData,panel=panel.xyplotWithDiag)

panel.xyplotWithLine <- function(x,y,intercept=0,slope=1,...) {
    panel.xyplot(x,y,...)
    panel.abline(a=intercept,b=slope,...)
}

xyplot(y~x, data=myData, panel=panel.xyplotWithLine)
xyplot(y~x, data=myData, 
    inter=0.5, slope=0, pch=16,
    lwd=2, col="gray30", lty=2,
    panel=panel.xyplotWithLine,
)
