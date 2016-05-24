plot.tgram <-
function(x, xlim=NULL, ylim=NULL, colores=c("red","green"),leyenda = c("lumen","double wall"), lwd=2, add=FALSE,traq.0=TRUE,bg.legend=NULL,...)
{

traq <- x$traq
puntos <- x$cut.points
que <- x$what

            s<- 1:(length(puntos[,1])-1)
if (add==FALSE) plot(traq, type="l", xlim=xlim, ylim=ylim,...) else lines(traq, ...)
if(traq.0 ==TRUE) lines(x$traq0, col="pink")
lines(traq,...)
points(puntos, col="red")
color <- colores[que]
segments(puntos[s,1],puntos[s,2],puntos[s+1,1],puntos[s+1,2], lwd=2,col=color )
legend(min(axis(1)),max(axis(2)),leyenda, col=colores, lwd=lwd,bg=bg.legend)
}

