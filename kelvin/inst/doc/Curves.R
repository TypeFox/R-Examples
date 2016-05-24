## ---- echo=TRUE----------------------------------------------------------
library(kelvin)

## ---- echo=TRUE, results='asis'------------------------------------------
plot_kelvin_curves <- function(){
  # Kelvin functions (order 0 by default)
  num.pts <- 1e3
  curve(expr=Kei,
        from=0.001, to=5, n=num.pts,
        ylim=c(-7,7), xlim=c(0,5), 
        yaxs="i",
        xaxs="i",
        main="Fundamental Kelvin functions",
        ylab="Ke(x) or Be(x)", lwd=2)
  curve(expr=Ker,from=0.001,to=5,n=num.pts,add=T,col='red', lwd=2)
  
  # complementary Kelvin functions (order 0 by default)
  curve(expr=Bei, from=0.001, to=5, n=num.pts, add=TRUE, lty=2, lwd=2)
  curve(expr=Ber, from=0.001, to=5 ,n=num.pts, add=TRUE, col='red', lty=2, lwd=2)
  legend(0.5, 5, c(expression(Kei[0]), expression(Ker[0])), col=c(1,2), lty=c(1,1), lwd=2)
  legend(2.8, -3.8, c(expression(Bei[0]), expression(Ber[0])), col=c(1,2), lty=c(2,2), lwd=2)
  
  xseq <- seq.int(0.001, 5, length.out=num.pts)
  
  Knu <- Keir(xseq, nSeq=6, return.list=FALSE)
  matplot(xseq, Re(Knu), type="l", xaxs="i", xlim=c(0,5), yaxs="i", ylim=c(-7,7),
          lty=1, lwd=2, main="Fundamental and higher order Kelvin functions (Ker)")
  legend(3.5, 7, 0:5, col=1:6, lty=1, lwd=2)
  
  Bnu <- Beir(xseq, nSeq=6, return.list=FALSE)
  matplot(xseq, Re(Bnu), type="l", xaxs="i", xlim=c(0,5), yaxs="i", ylim=c(-7,7),
          lty=1, lwd=2, main="Fundamental and higher order complimentary Kelvin functions (Ber)")
  legend(0.5, 7, 0:5, col=1:6, lty=1, lwd=2)
}

## ---- echo=TRUE, fig.show='hold', fig.width=8.5, fig.height=11-----------
plot_kelvin_curves()

