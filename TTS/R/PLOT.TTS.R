##' Time-Temperature Superposition (TTS) plots
##' 
##' Plots of TTS results: experimental data, horizontal and vertical shifts,
##' TTS data, TTS Master Curve fitting with B-Splines and bootstrap confidence
##' intervals are deployed.
##' 
##' TTS plots are performed from the outputs of TTS function: data, aT, bT,
##' TTS.data, TTS.gam y residuals.
##' 
##' @param x TTS object.
##' @return The following values are returned: \item{PLOT.data()}{Generic
##' function to plot the experimental data. By default log10.module versus
##' log10.frequency.} \item{PLOT.aT()}{Generic plot of the horizontal shifts
##' corresponding to each curve (modulus versus frequency) obtained on
##' temperature.} \item{PLOT.bT()}{Generic plot of the vertical shifts
##' corresponding to each curve (modulus versus frequency) obtained on
##' temperature.} \item{PLOT.TTS.data()}{Generic plot of the experimental data
##' horizontally and vertically shifted with respect to a the curve
##' corresponding to the reference temperature.} \item{PLOT.TTS.gam()}{Generic
##' plot of the Master Curve B-splines estimation with bootstrap confidence
##' intervals at 95 per cent.} \item{PLOT.res()}{Generic plot of the residuals
##' of Master Curve B-splines fitting.}
##' @author Antonio Meneses \email{antoniomenesesfreire@@hotmail.com}, Salvador
##' Naya \email{salva@@udc.es} and Javier Tarrio-Saavedra
##' \email{jtarrio@@udc.es}
##' @references Naya, S., Meneses A., Tarrio-Saavedra, J., Artiaga R.,
##' Lopez-Beceiro, J. and Gracia-Fernandez C. (2013) New method for estimating
##' shift factors in time-temperatura superposition models. Journal of Thermal
##' Analysis and Calorimetry. ISSN 1388-6150. DOI 10.1007/s10973-013-3193-1.\cr
##' 
##' Williams, M. L. (1964) Structural analysis of Viscoelastic materials. AIAA
##' Journal, 785-808.\cr
##' 
##' Artiaga R., Garcia A. Fundamentals of DMA. In: 'Thermal analysis.
##' Fundamentals and applications to material characterization' (ed.: Artiaga
##' R.) Publicaciones de la Universidade da Coruna, A Coruna, Spain, 183-206
##' (2005).\cr
##' @keywords TTS.PLOT()
##' @examples
##' 
## library(TTS)
##' ## TTS object applied to PC dataset.
##' data(PC)
##' Derive <- TTS(PC)
##' x <- Derive
##' ## Generic plots for TTS analysis
##' PLOT <- PLOT.TTS(x)
##' names(PLOT)
##' ##[1] "PLOT.data"     "PLOT.aT"       "PLOT.bT"       "PLOT.TTS.data"
##' ##[5] "PLOT.TTS.gam"  "PLOT.res"
##' ## Generic plots of: data, aT, bT, TTS.data, TTS.gam and res
##' PLOT$PLOT.data(main="PLOT: Data",xlab="log10.Frequency (rad/s)",ylab="log10.E'(Pa)")
##' PLOT$PLOT.aT(main="PLOT: horizontal translation factors", xlab="Temperature", ylab="aT")
##' PLOT$PLOT.bT(main="PLOT: vertical translation factors", xlab="Temperature",ylab="bT")
##' PLOT$PLOT.TTS.data(xlab="log10.Frequency (rad/s)",ylab="log10.E'(Pa)")
##' PLOT$PLOT.TTS.gam( xlab="log10.Frequency (rad/s)", ylab = "log10.E'(Pa)",
##' main = "Fitted gam, Bootstrap confidence intervals",
##' sub = "Reference temperature = 150 degrees celsius")
##' PLOT$PLOT.res(main="TTS: gam residual", xlab="Fitted", ylab="Standardized residuals")
##' 
##' @export PLOT.TTS
PLOT.TTS <-
function(x){
    if(x$data[1,2]<=x$data[length(x$data[,2][x$data[,3]==unique(x$data[,3])[1]]),2])
       { m <- 1
    }else m <- 2
  PLOT.data <- function(y=x$data,...){
      COL <- c(1:7,colors()[82:150])
      plot.data <- plot(subset(y[,c(1,2)],y[,3]==unique(y[,3])[1]),type="b",
                     xlim=c(min(y[,1]),max(y[,1])),ylim=c(min(y[,2]),max(y[,2])),
                     col=1,las=1,pch=20,cex=2,...)
      for(i in 2:length(unique(y[,3])))
      points(subset(y[,c(1,2)],y[,3]==unique(y[,3])[i]),type="b",col=COL[i],pch=20,cex=2)
      if(m==1) POS <- "bottomright" else POS <- "topright"
      legend(POS,c("Temperature",as.character(unique(y[,3]))),
      col=c(0,COL),lty=3,lwd=5,bty="n")
  }
  PLOT.aT <- function(y=x$aT,z=unique(x$data[,3]),...)
      plot(z,y,type="b", lwd=3,pch=20,las=1,col=4,...)
  PLOT.bT <- function(y=x$bT,z=unique(x$data[,3]),...)
      plot(z,y,type="b",lwd=3,pch=20,las=1,col=4,...)
  PLOT.TTS.data <- function(y=x$TTS.data,z=x$ref.temp,...){
      COL <- c(1:7,colors()[82:150])
      plot.data <- plot(subset(y[,c(1,2)],y[,3]==unique(y[,3])[1]),type="p",
                     xlim=c(min(y[,1]),max(y[,1])),ylim=c(min(y[,2]),max(y[,2])),
                     main=paste("PLOT: TTS\n Reference temperature = ",z),
                     col=1,las=1,pch=20,cex=2,...)
      for(i in 2:length(unique(y[,3])))
      points(subset(y[,c(1,2)],y[,3]==unique(y[,3])[i]),type="p",col=COL[i],pch=20,cex=2)
      if(m==1) POS <- "bottomright" else POS <- "topright"
      legend(POS,c("Temperature",as.character(unique(y[,3]))),
      col=c(0,COL),lty=3,lwd=5,bty="n")
  }
  PLOT.TTS.gam <- function(y=x$TTS.data,z=x$ref.temp,u=x$TTS.gam,v=x$I.lower,w=x$I.upper,...){
      COL <- c(1:7,colors()[82:150])
      plot.data <- plot(subset(y[,c(1,2)],y[,3]==unique(y[,3])[1]),type="p",
                     xlim=c(min(y[,1]),max(y[,1])),ylim=c(min(y[,2]),max(y[,2])),
                     col=1,las=1,pch=1,cex=1,lwd=3,...)
      for(i in 2:length(unique(y[,3])))
      points(subset(y[,c(1,2)],y[,3]==unique(y[,3])[i]),type="p",col=COL[i],pch=1,cex=1,lwd=3)
      if(m==1) POS <- "bottomright" else POS <- "topright"
      legend("bottomright",c("Temperature",as.character(unique(y[,3]))),
      col=c(0,COL),lty=3,lwd=5,bty="n")
      points(u,col=4,type="l",lwd=3)
      lines(u[,1],v,lty=2,col=2)
      lines(u[,1],w,lty=2,col=2)
  }
  res.stand <- (x$residuals-mean(x$residuals))/sd(x$residuals)
  PLOT.res <- function(y=res.stand,z=x$TTS.gam,...){
      plot(z[,2],y,type="p", cex=1.5, pch=20, las=1, col=4,...)
      abline(h=0,lwd=2,col=2)
  }

list(PLOT.data=PLOT.data, PLOT.aT=PLOT.aT, PLOT.bT=PLOT.bT, PLOT.TTS.data=PLOT.TTS.data,
     PLOT.TTS.gam=PLOT.TTS.gam,PLOT.res=PLOT.res)
}


