plot.MaxEntICA.BinBin <- function(x, ICA.Fit, Type="Density", Xlab, col, Main, ...){
Object <- x 

if (missing(Xlab)) {Xlab <- expression(R[H]^2)}
if (missing(col)) {col <- c(8)}
if (missing(Main)) {Main <- " "}  
    
    if (Type=="Density"){    

    plot(density(ICA.Fit$R2_H, na.rm = T), xlab=Xlab, ylab="Density", main=Main, lwd=2, col=col, ...)
    abline(v=Object$R2_H, lwd=2) }

    if (Type=="Freq"){
      
    h <- hist(ICA.Fit$R2_H, plot = FALSE)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=ICA.Fit$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
        plot(h,freq=T, xlab=Xlab, ylab="Frequency", col=col, main=Main, ...)
        abline(v=Object$R2_H, lwd=2) }
}