### R code from vignette source 'Fig11GriffisStedinger.Rnw'

###################################################
### code chunk number 1: Fig11GriffisStedinger.Rnw:44-74
###################################################
tau4_LP3 <- function (tau3_LP3) {
 # the approximations yield values of tau4 accurate to 
 #  within 0.008 over the range tau3min <= tau3 <= 0.9
 gammax <- c(-1.4, -1, -0.5, 0, 0.5, 1, 1.4)
 table3abcd <- matrix(c(0.0602, -0.1673,  0.8010,  0.2897,
                        0.0908, -0.1267,  0.7636,  0.2562,
                        0.1166, -0.0439,  0.6247,  0.2939,
                        0.1220,  0.0238,  0.6677,  0.1677,
                        0.1152,  0.0639,  0.7486,  0.0645,
                        0.1037,  0.0438,  0.9327, -0.0951,
                        0.0776,  0.0762,  0.9771, -0.1394),
                        ncol=4, nrow=7, byrow=TRUE,
                        dimnames=list(paste("gammax=", gammax, sep=""), c("a","b","c","d")))

 tau3s <- matrix(c(rep(1, length(tau3_LP3)), tau3_LP3, tau3_LP3^2, tau3_LP3^3), 
                 nrow=4, ncol=length(tau3_LP3), byrow=TRUE,
                 dimnames=list(c("1","tau3","tau3^2","tau3^3"), 
                            paste("tau3=", tau3_LP3, sep="")))

 tau3min <- c(-0.2308, -0.1643, -0.0740, 0, 0.0774, 0.1701, 0.2366) %*% 
            t(rep(1, length(tau3_LP3)))
 tau3max <- rep(0.9, 7) %*% t(rep(1, length(tau3_LP3)))
 tau3s2 <- rep(1,7) %*% t(tau3_LP3)
 keeptau4 <- (tau3s2 >= tau3min*0.9999) & (tau3s2 <= tau3max*1.0001)

 tau4s <- table3abcd %*% tau3s
 tau4s[!keeptau4] <- NA

 return(tau4s)
}


###################################################
### code chunk number 2: Fig11GriffisStedinger.Rnw:78-89
###################################################
gammax <- c(-1.4, -1, -0.5, 0, 0.5, 1, 1.4)
tau3 <- seq(-0.2, 1, by=0.1)
tau4 <- tau4_LP3(tau3)
plot(c(-0.4, 1), c(0, 1), type="n", xlab="L-skewness, tau3", ylab="L-Kurtosis, tau4")
grid()
for (i in 1:7) {
 lines(tau3, tau4[i,], type="b", lty=i, pch=i)
}
legend("topleft", legend=c(gammax, "OLB"), title=expression(gamma[x]), 
       pch=c(1:7,NA), lty=c(1:7,1), lwd=c(rep(1,7),2))
curve((5*x^2 - 1)/4, add=TRUE, lwd=2)


###################################################
### code chunk number 3: Fig11GriffisStedinger.Rnw:97-98
###################################################
library(nsRFA)


###################################################
### code chunk number 4: Fig11GriffisStedinger.Rnw:101-139
###################################################
t3b=-0.2; t3t=0.9; t4b=-0.1; t4t=0.8
 plot(c(t3b+0.05, t3t-0.05), c(t4b, t4t), type="n", 
      xlab=expression(tau[3]), ylab=expression(tau[4]), main="")
  grid()
tipi <- c(1,2,1,4,5,6)
spessori <- c(1,1.3,1.1,1.3,1.1,1.1)
colori <- c(1,1,"darkgrey",1,1,1)
tau3 <- seq(0, 0.9, by=0.1)
tau4 <- tau4_LP3(tau3)
tau4r <- apply(tau4, 2, range, na.rm=TRUE)
 polygon(x=c(tau3, rev(tau3)), y=c(tau4r[1,], rev(tau4r[2,])), 
         density=20, col="darkgrey", border="darkgrey", angle=-45)

 GPA <- function(x) 0.20196*x + 0.95924*x^2 - 0.20096*x^3 + 0.04061*x^4
  curve(GPA, t3b, t3t, add=TRUE, lty=tipi[5], lwd=spessori[5])
 GEV <- function(x) {
   0.10701 + 0.1109*x + 0.84838*x^2 - 0.06669*x^3 + 
    0.00567*x^4 - 0.04208*x^5 + 0.03763*x^6
  }
  curve(GEV, t3b, t3t, add=TRUE, lty=tipi[4], lwd=spessori[4])
 GLO <- function(x) 0.16667 + 0.83333*x^2
  curve(GLO, t3b, t3t, add=TRUE, lty=tipi[6], lwd=spessori[6])
 LN3 <- function(x) 0.12282 + 0.77518*x^2 + 0.12279*x^4 - 0.13638*x^6 + 0.11368*x^8
  curve(LN3, t3b, t3t, add=TRUE, lty=tipi[1], lwd=spessori[1])
 PE3 <- function(x) 0.1224 + 0.30115*x^2 + 0.95812*x^4 - 0.57488*x^6 + 0.19383*x^8
  curve(PE3, t3b, t3t, add=TRUE, lty=tipi[2], lwd=spessori[2])

  points(0, 0, pch=3, cex=1.2)
  points(0, 0.1226, pch=2, cex=1.2)
  points(1/3, 1/6, pch=5, cex=1.2)
  points(0.1699, 0.1504, pch=6, cex=1.2)
  points(0, 1/6, pch=4, cex=1.2)
  curve((5*x^2 - 1)/4, t3b, t3t, add=TRUE, lwd=2)

  legend("bottomright", c("EXP", "EV1", "LOG", "NOR", "UNIF"), 
         pch=c(5, 6, 4, 2, 3), bty="n")
  legend("topleft", legend=c("LN3","P3","LP3","GEV","GP","GL","OLB"), 
         lty=c(tipi,1), lwd=c(spessori,2), col=c(colori,1), bty="n")


