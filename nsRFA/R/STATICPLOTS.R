Lmoment.ratio.diagram <- function(grid=TRUE, ...) {

  # INPUT
  # ...
  # OUTPUT
  # plot dell'Lmoment ratio diagram con le curve corrispondenti a:
  # E = exponential
  # G = gumbel
  # L = logistic
  # N = normal
  # U = uniform
  # GLO = generalized logistic
  # GEV = generalized extreme-value
  # GPA = generalized Pareto
  # LN3 = lognormal
  # PE3 = Pearson type III

  plot(c(-0.2,0.6),c(-0.1,0.4),type="n",xlab="L-CA",ylab="L-kur",main="", ...)
  if (grid) {
   abline(v=seq(-0.2,0.6,by=0.02),col="lightgray",lty = "dotted")
   abline(h=seq(-0.1,0.4,by=0.02),col="lightgray",lty = "dotted")
  }

  GPA <- function(x) 0.20196*x + 0.95924*x^2 - 0.20096*x^3 + 0.04061*x^4
  curve(GPA, -0.2, 0.6, add=TRUE, lty=2)

  GEV <- function(x) 0.10701 + 0.11090*x + 0.84838*x^2 - 0.06669*x^3 + 0.00567*x^4 - 0.04208*x^5 + 0.03763*x^6
  curve(GEV, -0.2, 0.6, add=TRUE, lty=5)

  GLO <- function(x) 0.16667 + 0.83333*x^2
  curve(GLO, -0.2, 0.6, add=TRUE, lty=3)

  LN3 <- function(x) 0.12282 + 0.77518*x^2 + 0.12279*x^4 - 0.13638*x^6 + 0.11368*x^8
  curve(LN3, -0.2, 0.6, add=TRUE, lty=4)

  PE3 <- function(x) 0.12240 + 0.30115*x^2 + 0.95812*x^4 - 0.57488*x^6 + 0.19383*x^8
  curve(PE3, -0.2, 0.6, add=TRUE, lty=1)

  #U <- c(0,0)
  points(0,0, pch=3, cex=1.2)

  #N <- c(0,0.1226)
  points(0,0.1226, pch=2, cex=1.2)

  #E <- c(1/3,1/6)
  points(1/3,1/6, pch=5, cex=1.2)

  #G <- c(0.1699,0.1504)
  points(0.1699,0.1504, pch=6, cex=1.2)

  #L <- c(0,1/6)
  points(0,1/6, pch=4, cex=1.2)

  legend("bottomright", c("E","G","L","N","U"),pch=c(5,6,4,2,3),xjust=1,yjust=0,bty="n")
  legend("topleft", c("GLO","GEV","GPA","LN3","PE3"),lty=c(3,5,2,4,1),xjust=0,yjust=1,bty="n")
}


# -------------------------------------------------------------------------------------------- #

Lspace.HWvsAD <- function (grid=TRUE, ...) {

  plot(c(-0.1,0.5),c(0.1,0.6),type="n",xlab="L-CA",ylab="L-CV", ...)
  if(grid) {grid()}

  abline(-0.2,1,lty=2)
  abline(0.4,1,lty=2)

  lines(c(0.23,0.23),c(0.08,0.62),lty=1,lwd=2)

  text(0.05,0.23,expression(HW),cex=1.5)
  text(0.37,0.4,expression(AD),cex=1.5)

}


# ------------------------------------------------------------------------------------------ #

Lspace.limits <- function (grid=TRUE, ...) {

  # Produce una figura dello spazio di interesse

  plot(c(-0.2,0.8),c(0.1,0.8),type="n",xlab="L-CA",ylab="L-CV", ...)
  if(grid) {grid()}

  abline(0.5,0.5,lty=3, col=1)
  abline(-0.2,1,lty=2, col=1)
  abline(0.4,1,lty=2, col=1)

  lines(c(-0.1,0.3,0.5,0.5,0.2,-0.1,-0.1),c(0.1,0.1,0.3,0.6,0.6,0.3,0.1))

  legend(0.8,0.1,c("a","b","c"),lty=c(3,2,1),col=c(1,1,1),bty="n",xjust=1,yjust=0)
}

