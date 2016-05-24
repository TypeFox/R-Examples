barchart <- function(x,col,my.title,pareto=FALSE,freq.cumul=FALSE,family="Courier") {
# Save default values for par()
  sauve.par <- par(no.readonly = TRUE)
#par(tck=1)
  tmp <- par()$mai
  tmp[1] <- par()$pin[2]/5+0.1
  par(mai=tmp)
  if (missing(my.title)) {
    pareto.my.title <- if (pareto) "Pareto" else "Bar"
    my.title <- paste(paste(pareto.my.title,"chart"),paste("of variable",
                                                     deparse(substitute(x))),sep="\n")
    cumul.my.title <- "and cumulative frequencies"
    if (freq.cumul) {
      my.title <- paste(paste(pareto.my.title,"chart"),cumul.my.title,
                     paste("of variable",deparse(substitute(x))),sep="\n")
    }
  }
  n <- length(x)
  lx <- levels(x)
  nx <- length(lx)
  rangs <- if (pareto) order(table(x),decreasing=TRUE) else 1:nx
  tx <- table(x)[rangs]
  if (freq.cumul) tx <- tx/n
  if (freq.cumul) ylim <- c(0,1) else ylim <- c(0,max(tx))
  posy <- if (freq.cumul) -0.1 else -max(tx)/10
  if (missing(col)) col <- 1:nx
  col <- col[rangs]
  spaces <- c((1-((nx-1)*0.1/10+nx*0.1))/2*1000,rep(10*1/nx*20,nx-1))/100
  spaces2 <- spaces
  spaces2[1] <- spaces2[1]+0.12
  if (nx==2) spaces <- rev(spaces)
  if (nx==2) spaces2 <- rev(spaces2)
  plot.new()
  plot.window(xlim=c(0,1),ylim=ylim)
  
  for (i in seq(from=ylim[2]/5,to=ylim[2],by=ylim[2]/5)) {
    polygon(c(-0.04,1.1,1.1,-0.04),c(i-ylim[2]/20,i-ylim[2]/20,i,i),col=rgb(211/256,211/256,211/256,0.3),border=NA)
    abline(h=c(i-ylim[2]/20,i),col="lightgrey")
  }
  
  barplot(tx,col="darkgray",width=rep(0.1,n),space=spaces2,xlim=c(0,1),
          border="white",axis.lty=0,names.arg="",axes=FALSE,ylim=ylim,add=TRUE)
  axis(2,col.ticks="gray")
  r <- barplot(tx,col=col,width=rep(0.1,n),space=spaces,xlim=c(0,1),
               border="black",axis.lty=0,names.arg="",
               add=TRUE,axes=FALSE)
  if (freq.cumul) {
    eps <- 0.02
    points(r,cumsum(tx),type="l",col="#BF0000",lwd=3,lend="square") # red line
    points(r+eps,cumsum(tx)-eps,type="l",col=rgb(0,0,0,0.2),lwd=2,ljoin="bevel") # black shadow
    points(r,cumsum(tx),type="p",pch=20,col="#BF0000") # red disc
  }
  
  title(my.title,family=family)
  arrowaxis(x=FALSE,y=TRUE) 
  extr <- par("usr")
  if (freq.cumul)  arrows(0,0,extr[2]+0.05*diff(extr[1:2]),0,xpd=TRUE,length=0.1)
  text(r,posy,lx[rangs],srt=45,xpd=TRUE,family="Hershey",vfont=c("script","plain"),cex=1.2)
  abline(h=0)
  par(sauve.par)
}
