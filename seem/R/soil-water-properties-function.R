soil.water <- function(theta, param, title="", plot=F, pdfout=F) {

# Based on Brooks Corey 
# arguments:
# theta:  seq of soil water content as wetness (vol fraction)
# param:
#  Pb bubbling suction in cm
#  b pore size distribution index (adim) 
#  theta.r residual water
#  theta.s total porosity
# option T or F to plot

Pb <- param[1]; b <- param[2]
theta.r <- param[3]; theta.s <- param[4]

# effective saturation
theta.e <- (theta-theta.r)/(theta.s-theta.r)

# retention curve
# P is matric suction in cm
P <- Pb/(theta.e)^b

# capacity
C <- (1/b)*(Pb/P)^(1/b-1)*(-Pb/P^2)

# Hydraulic conductivity in cm/sec
# saturated conductivity
Ks <- 86/((b+1)*(2*b+1))*((theta.s-theta.r)/Pb)^2
# Green-Ampt wetting front suction
Pf <- ((2*b+3)/(b+3))*Pb

# non saturated
K <- Ks*theta.e^(2*b+3)

# logarithmic
LP <- log(P)
LC <- log(-C)
LK <- log(K)

if (plot==T) {
mat<- matrix(1:4,2,2,byrow=T)
layout(mat,c(3.5,3.5),c(3.5,3.5),respect=TRUE)
par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

 xlabel= "Theta content"
 plot(theta, theta.e, type="l", xlim=c(0,1.3*max(theta)), ylim=c(0,1),ylab="Theta-e",xlab=xlabel)
 title(title, cex.main=0.7) 
 plot(theta, P,  type="l", xlim=c(0,1.3*max(theta)), ylim=c(0,max(P)), ylab="P [cm]",xlab=xlabel)
 title(title, cex.main=0.7) 
 plot(theta, C,  type="l", xlim=c(0,1.3*max(theta)), ylim=c(min(C),max(C)), ylab="C [1/cm]",xlab=xlabel)
 title(title, cex.main=0.7) 
 plot(theta, K,  type="l", xlim=c(0,1.3*max(theta)), ylim=c(0,max(K)), ylab="K [cm/s]",xlab=xlabel)
 title(title, cex.main=0.7) 

 if(pdfout==F) {
  win.graph()
  mat<- matrix(1:4,2,2,byrow=T)
  layout(mat,c(3.5,3.5),c(3.5,3.5),respect=TRUE)
  par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
 }
 plot(theta, theta.e, type="l", xlim=c(0,1.3*max(theta)), ylim=c(0,1),ylab="Theta-e",xlab=xlabel)
 title(title, cex.main=0.7) 
 par(ylog=T)
 plot(theta, P,  type="l", log="y",xlim=c(0,1.3*max(theta)), ylab="P [cm] log scale",xlab=xlabel)
 title(title, cex.main=0.7) 
 plot(theta, -C, type="l", log="y",xlim=c(0,1.3*max(theta)), ylab="C [1/cm] log scale",xlab=xlabel)
 title(title, cex.main=0.7) 
 plot(theta, K,  type="l", log="y",xlim=c(0,1.3*max(theta)), ylab="K [cm/s] log scale",xlab=xlabel)
 title(title, cex.main=0.7) 

}
theta.e <- round(theta.e,2)
P <- round(P)
C <- signif(C,2)
Pf <- round(Pf)
Ks <- signif(Ks,2)
K <- signif(K,2)
LP <- round(LP)
LC <- round(LC)
LK <- round(LK)

out <- list(Ks=Ks, Pf=Pf, var.theta=data.frame(theta,theta.e,P,C,K,LP,LC,LK))
return(out)
}
