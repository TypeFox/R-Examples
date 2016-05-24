plot.GoodPretreatContCont <- function(x, main, col, Type="Percent", Labels=FALSE, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), ...) {
  Object <- x 
  if (missing(main)) {main = " "}
  if (missing(col)) {col=8}
  
  dat <- data.frame(cbind(Object$T0T1, Object$Rho2.Min))
  colnames(dat) <- c("T0T1", "Rho2.Min")
  dev.new()
  par=Par
  
  if (Type=="Freq"){
    
    h <- hist(dat$Rho2.Min, ...)
    h$density <- h$counts/sum(h$counts)
    cumulMidPoint <- ecdf(x=dat$Rho2.Min)(h$mids)
    labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
    
    if (Labels==FALSE){
      plot(h, freq=T, xlab=expression(rho[min]^2), ylab="Frequency", col=col, main=main)
    }
    if (Labels==TRUE){
      plot(h, freq=T, xlab=expression(rho[min]^2), ylab="Frequency", col=col, main=main, labels=labs)}
  }
  
  if (Type=="Percent"){
    
    h <- hist(dat$Rho2.Min, ...)
    h$density <- h$counts/sum(h$counts)
    cumulMidPoint <- ecdf(x=dat$Rho2.Min)(h$mids)
    labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
    
    if (Labels==FALSE){
      plot(h, freq=F, xlab=expression(rho[min]^2), ylab="Percentage", col=col, main=main)}
    if (Labels==TRUE){
      plot(h, freq=F, xlab=expression(rho[min]^2), ylab="Percentage", col=col, main=main, labels=labs)}
  }
  
  if (Type=="CumPerc"){
    h <- hist(dat$Rho2.Min, breaks=length(dat$Rho2.Min), ...)
    h$density <- h$counts/sum(h$counts)
    cumulative <- cumsum(h$density)
    plot(x=h$mids, y=cumulative, xlab=expression(rho[min]^2), ylab="Cumulative percentage", col=0, main=main)
    lines(x=h$mids, y=cumulative)
  }
  
  
}
