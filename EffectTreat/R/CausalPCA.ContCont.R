CausalPCA.ContCont <- function(x, Min=-1, Max=1, Cex.Letters=3, Cex.Corrs=2, Lines.Rel.Width=TRUE, Col.Pos.Neg=TRUE) {

if (is(x, "PCA.ContCont")==TRUE){  

dat <- cbind(x$Pos.Def, x$PCA)
colnames(dat) <- c("T0T1", "T0S", "T1S", "PCA")
sub <- dat[dat$PCA >= Min & dat$PCA <= Max,] 
 
med_T0T1 <- round(median(sub$T0T1), digits=2)
med_T0S <- unique(round(sub$T0S, digits=2))
med_T1S <- unique(round(sub$T1S, digits=2))

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(0:10, 0:10, axes=F, xlab="", ylab="", type="n")  
par(oma=c(0, 0, 0, 0))
text(1, 1, expression(S), cex=Cex.Letters)
text(9, 9, expression(T[0]), cex=Cex.Letters)
text(9, 1, expression(T[1]), cex=Cex.Letters)
  
text(5, 6, med_T0S, cex=Cex.Corrs)
text(5, 0.5, med_T1S, cex=Cex.Corrs)
text(9.6, 5, med_T0T1, cex=Cex.Corrs)

  
  if (Lines.Rel.Width==TRUE){
  
    if (Col.Pos.Neg==FALSE) {col_T1S <- col_T0S <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_T1S <- col_T0S <- col_T0T1 <- 1
      if (med_T1S<0) {col_T1S <- "red"}
      if (med_T0S<0) {col_T0S <- "red"}
      if (med_T0T1<0) {col_T0T1 <- "red"}
    }
  
  segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1+(abs(med_T0S)*5), col=col_T0S)
  segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1+(abs(med_T1S)*5), col=col_T1S)
  segments(x0=9, y0=8, x1=9, y1=2, lwd=1+(abs(med_T0T1)*5), col=col_T0T1)
  
  }

  if (Lines.Rel.Width==FALSE){
    
    if (Col.Pos.Neg==FALSE) {col_T1S <- col_T0S <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_T1S <- col_T0S <- col_T0T1 <- 1
      if (med_T1S<0) {col_T1S <- "red"}
      if (med_T0S<0) {col_T0S <- "red"}
      if (med_T0T1<0) {col_T0T1 <- "red"}
    }
    
  segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1, col=col_T0S)
  segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1, col=col_T1S)
  segments(x0=9, y0=8, x1=9, y1=2, lwd=1, col=col_T0T1)
  }

}
  

if (is(x, "Multivar.PCA.ContCont")==TRUE){  
  stop("x should be fitted object of class PCA.ContCont")
 }
}  
  


