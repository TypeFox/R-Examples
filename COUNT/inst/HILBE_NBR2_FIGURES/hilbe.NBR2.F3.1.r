# hilbe_NBR2_F3_1.r
# Table 3.1
# Script for creating Poisson and negative binomial (NB2) distributions 
#  on single graph with user specified number of observations (obs), 
#  mean (mu), counts (0:specified), and NB2 ancillary parameter value (alpha).
#  Additional enhanced graphic can be displayed; -ggplot2- must be installed. 
#  From Hilbe, Negative Binomial Regression, 2nd ed, Cambridge Univ. Press
#  Table 3.1; Figure 3.1
#
# ------------------------------------------
# User specified values; defaults displayed
obs   <- 11
mu    <- 2
y     <- 0:10
alpha <- 1.5
# ------------------------------------------
amu <- mu*alpha
ynb2 = exp(
        y*log(amu/(1+amu)) 
     - (1/alpha)*log(1+amu) 
     + log( gamma(y +1/alpha) )
     - log( gamma(y+1) ) 
     - log( gamma(1/alpha) ) 
)
xbar = "mu"
a = "alpha"
plot(  y, ynb2, col="red", pch=5,
  main="Poisson vs Negative Binomial PDFs")
  lines( y, ynb2, col="red")
  points(y, yp2,  col="blue", pch=2)
  lines( y,  yp2, col="blue")
  legend(4.3,.40, 
    c("Negative Binomial: mean=2, a=1.5",
      "Poisson: mean=2"),
    col=( c("red","blue") ),
    pch=( c(5,2) ),
    lty=1)
#========= FOR NICER GRAPHIC =======================
zt <- 0:10           #zt is Zero to Ten
x <- c(zt,zt)        #two zt's stacked for use with ggplot2
newY <- c(yp2, ynb2) #Now stacking these two vars
Distribution <- gl(n=2,k=11,length=22,
  label=c("Poisson","Negative Binomial")
)
NBPlines <- data.frame(x,newY,Distribution)
library("ggplot2")
ggplot( NBPlines, aes(x,newY,shape=Distribution,col=Distribution ) ) +
  geom_line() + geom_point()
