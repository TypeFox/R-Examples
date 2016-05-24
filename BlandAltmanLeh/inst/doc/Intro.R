## ------------------------------------------------------------------------
A <- c(-0.358, 0.788, 1.23, -0.338, -0.789, -0.255, 0.645, 0.506, 
0.774, -0.511, -0.517, -0.391, 0.681, -2.037, 2.019, -0.447, 
0.122, -0.412, 1.273, -2.165)
B <- c(0.121, 1.322, 1.929, -0.339, -0.515, -0.029, 1.322, 0.951, 
0.799, -0.306, -0.158, 0.144, 1.132, -0.675, 2.534, -0.398, 0.537, 
0.173, 1.508, -1.955)

## ---- cache=FALSE, results='hide'----------------------------------------
plot(A, B, main="Scatter plot")
abline(0,1)

## ---- results='hide'-----------------------------------------------------
plot((A+B)/2, A-B, main="Mean-Difference-Plot")

## ---- results='hide'-----------------------------------------------------
library(BlandAltmanLeh)
bland.altman.plot(A, B, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences")

## ------------------------------------------------------------------------
library(ggplot2)
pl <- bland.altman.plot(A, B, graph.sys = "ggplot2")
print(pl)

## ---- results='hide'-----------------------------------------------------
A <- rnorm(50)
B <- A +runif(50, -.3, .3)
pl <- bland.altman.plot(A, B, graph.sys="ggplot2", conf.int=.95)
print(pl)

# or in base-graphics:
bland.altman.plot(A, B, conf.int=.95, pch=19)

## ---- results='hide'-----------------------------------------------------
A <- c(-0.358, 0.788, 1.23, -0.338, -0.789, -0.255, 0.645, 0.506, 
0.774, -0.511, -0.517, -0.391, 0.681, -2.037, 2.019, -0.447, 
0.122, -0.412, 1.273, -2.165)
B <- c(0.121, 1.322, 1.929, -0.339, -0.515, -0.029, 1.322, 0.951, 
0.799, -0.306, -0.158, 0.144, 1.132, -0.675, 2.534, -0.398, 0.537, 
0.173, 1.508, -1.955)
sex <- c( 1,1,1,1,2,2,2,1,1,1,2,2,2,2,2,1,1,2,1,2)

ba.stats <- bland.altman.stats(A, B)

plot(ba.stats$means, ba.stats$diffs, col=sex, 
     sub=paste("critical difference is", round(ba.stats$critical.diff,4)),
     main="make your own graph easily", ylim=c(-1.5,1.5), pch=18-sex)
abline(h = ba.stats$lines, lty=c(2,3,2), col=c("lightblue","blue","lightblue"), 
       lwd=c(3,2,3))
legend(x = "topright", legend = c("male","female"), fill = 1:2)

## ---- results='hide'-----------------------------------------------------
A <- c(7, 8, 4, 6, 4, 5, 9, 7, 5, 8, 1, 4, 5, 7, 3, 4, 4, 9, 3, 3, 
1, 4, 5, 6, 4, 7, 4, 7, 7, 5, 4, 6, 3, 4, 6, 4, 7, 4, 6, 5, 1, 1, 1, 1, 1, 1)
B <- c(8, 7, 4, 6, 3, 6, 9, 8, 4, 9, 0, 5, 5, 9, 3, 5, 5, 8, 3, 3, 
1, 4, 4, 7, 4, 8, 3, 7, 7, 5, 6, 7, 3, 3, 7, 3, 6, 5, 9, 5, 1, 1, 1, 1, 1, 1)

bland.altman.plot(A, B)

## ---- results='hide'-----------------------------------------------------
bland.altman.plot(A, B, sunflower=TRUE)

## ---- results='hide'-----------------------------------------------------
print( bland.altman.plot(A, B, graph.sys = "ggplot2", geom_count = TRUE) )

## ---- echo=FALSE---------------------------------------------------------
a<- rnorm(150)+0.6
b<- .02*a+.3*rnorm(150)

## ------------------------------------------------------------------------
library(ggExtra)
print(ggMarginal(bland.altman.plot(a, b, graph.sys = "ggplot2"),
           type = "histogram", size=4))

