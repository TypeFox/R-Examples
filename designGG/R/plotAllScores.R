# Function name: plotAllScores
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


plotAllScores <- function(plot.obj,fileName=NULL)
{
    #This is used in designGG main function
    scores            <-  plot.obj$scores
    cooling           <-  plot.obj$cooling
    startTemp         <-  plot.obj$startTemp
    temperature       <-  plot.obj$temperature
    temperature.step  <-  plot.obj$temperature.step
    nIterations       <-  plot.obj$nIterations
    optimality        <-  plot.obj$optimality
    sciNotation <- function(x, digits = 1) {
      if (!x) return(0)
      exponent <- floor(log10(x))
      base <- round(x / 10^exponent, digits)     
      return(as.expression(substitute(base * 10^exponent, 
        list(base = base, exponent = exponent)))) 
        }
    if(!is.null(fileName)){                           
    png(filename=paste(fileName,"SAplot.png", sep=""), width=580, height=480,
        bg=" light gray")}

    par( mfrow=c(2,1) )
    par( mai=c(0.2,0.95, 0.76, 0.39) )
    ylim <- range(scores,finite=TRUE,na.rm=TRUE)

    plot( 1:length(scores), scores, type="l", lwd=2, xaxt="n", xlab=NULL,
            ylab="Socres", col="blue", main="Simulated Annealing" )
            
    text( x=0.8*length(scores),
          y=ylim[1]+0.95*(ylim[2]-ylim[1]),
          paste("T0=", startTemp , ", T.end=", sciNotation( temperature ), sep="") )

    text( x=0.9*length(scores),
          y=ylim[1]+0.85*(ylim[2]-ylim[1]),
          paste("T step=", round(temperature.step,digits=3), sep="") )

    text( x=0.9*length(scores),
          y=ylim[1]+0.75*(ylim[2]-ylim[1]),
          paste("nIterations=", nIterations, sep="") )

    par( mai=c(0.95, 0.95, 0.2, 0.39) )

    plot( 1:length(cooling),col="blue", cooling, type="l", lwd=2,
            xlab="Times of Moving", ylab="Cooling")
    if(!is.null(fileName))    dev.off()
}
