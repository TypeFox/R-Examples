CausalDiagramContCont <- function(x, Min=-1, Max=1, Cex.Letters=3, Cex.Corrs=2, Lines.Rel.Width=TRUE, Col.Pos.Neg=TRUE) {
  
if (class(x)=="ICA.ContCont"){

  dat <- cbind(x$Pos.Def, x$ICA)
  colnames(dat) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "ICA")
  sub <- dat[dat$ICA >= Min & dat$ICA <= Max,] 
}
 
if (class(x)=="MICA.ContCont"){
  dat <- cbind(x$Pos.Def, x$MICA)
  colnames(dat) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "MICA")
  sub <- dat[dat$MICA >= Min & dat$MICA <= Max,] 
}  

  med_T0T1 <- round(median(sub$T0T1), digits=2)
  med_T0S0 <- unique(round(mean(sub$T0S0, na.rm=TRUE), digits=2))
  med_T0S1 <- round(median(sub$T0S1), digits=2)
  med_T1S0 <- round(median(sub$T1S0), digits=2)
  med_T1S1 <- unique(round(mean(sub$T1S1, na.rm=TRUE), digits=2))
  med_S0S1 <- round(median(sub$S0S1), digits=2)
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot(0:10, 0:10, axes=F, xlab="", ylab="", type="n")  
  par(oma=c(0, 0, 0, 0))
  text(1, 9, expression(S[0]), cex=Cex.Letters)
  text(1, 1, expression(S[1]), cex=Cex.Letters)
  text(0.4, 5, med_S0S1, cex=Cex.Corrs)
  text(9, 9, expression(T[0]), cex=Cex.Letters)
  text(9, 1, expression(T[1]), cex=Cex.Letters)
  text(4, 6.9, med_T1S0, cex=Cex.Corrs)
  text(4, 3.1, med_T0S1, cex=Cex.Corrs)
  text(5, 9.5, med_T0S0, cex=Cex.Corrs)
  text(5, 0.5, med_T1S1, cex=Cex.Corrs)
  text(9.6, 5, med_S0S1, cex=Cex.Corrs)
  
  
  if (Lines.Rel.Width==TRUE){
  
    if (Col.Pos.Neg==FALSE) {col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1
      if (med_S0S1<0) {col_S0S1 <- "red"} 
      if (med_T1S0<0) {col_T1S0 <- "red"}
      if (med_T0S1<0) {col_T0S1 <- "red"}
      if (med_T0S0<0) {col_T0S0 <- "red"} 
      if (med_T1S1<0) {col_T1S1 <- "red"} 
      if (med_T0T1<0) {col_T0T1 <- "red"}
    }
  
  segments(x0=1, y0=8, x1=1, y1=2, lwd=1+(abs(med_S0S1)*5), col=col_S0S1)
  segments(x0=1.5, y0=8, x1=8.5, y1=2, lwd=1+(abs(med_T1S0)*5), col=col_T1S0)
  segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1+(abs(med_T0S1)*5), col=col_T0S1)
  segments(x0=1.5, y0=9, x1=8.5, y1=9, lwd=1+(abs(med_T0S0)*5), col=col_T0S0)
  segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1+(abs(med_T1S1)*5), col=col_T1S1)
  segments(x0=9, y0=8, x1=9, y1=2, lwd=1+(abs(med_T0T1)*5), col=col_T0T1)
  
  }

  if (Lines.Rel.Width==FALSE){
    
    if (Col.Pos.Neg==FALSE) {col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1
      if (med_S0S1<0) {col_S0S1 <- "red"} 
      if (med_T1S0<0) {col_T1S0 <- "red"}
      if (med_T0S1<0) {col_T0S1 <- "red"}
      if (med_T0S0<0) {col_T0S0 <- "red"} 
      if (med_T1S1<0) {col_T1S1 <- "red"} 
      if (med_T0T1<0) {col_T0T1 <- "red"}
    }
    
    segments(x0=1, y0=8, x1=1, y1=2, lwd=1, col=col_S0S1)
    segments(x0=1.5, y0=8, x1=8.5, y1=2, lwd=1, col=col_T1S0)
    segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1, col=col_T0S1)
    segments(x0=1.5, y0=9, x1=8.5, y1=9, lwd=1, col=col_T0S0)
    segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1, col=col_T1S1)
    segments(x0=9, y0=8, x1=9, y1=2, lwd=1, col=col_T0T1)
  }

}  
  


