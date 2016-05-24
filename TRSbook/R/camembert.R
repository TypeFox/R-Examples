camembert <- function(x,col=NULL,family="Courier") {
  var <- sort(table(x)/length(x))
  sauve.par <- par(no.readonly = TRUE)
 #par(bg="#4B5475",fg="white")
  par(bg="white",fg="black")
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1] > pin[2])
    xlim <- (pin[1]/pin[2]) * xlim
  else ylim <- (pin[2]/pin[1]) * ylim
  
  plot.window(xlim, ylim, "", asp = 1)
  abline(h=(-2:2)*0.4,col="#7A86AF",lty=8)
  par(new=TRUE)
  pie(var,labels=paste(names(var)," (",format(100*round(var,2))," %)",sep=""),
      col=col,cex=0.7,border="black")
  title(main = paste("Pie chart",paste("for variable"
          ,deparse(substitute(x))),sep="\n"),family=family,cex=0.5,col.main="black")
 #par(sauve.par)
  box("outer", col = "white", lwd=7)
 #box("inner", col = "grey", lwd=2)
  t <- seq(0,1,length=100)
  t2p <- -2*pi * t + 90 * pi/180
  x <- c(0.82 * cos(t2p),rev(0.8 * cos(t2p)))
  y <- c(0.82 * sin(t2p),rev(0.8 * sin(t2p)))
  polygon(x, y, border = rgb(0,0,0,alpha=0) , col = rgb(0,0,0,alpha=0.3))
}
