skeletonsArchetypal <- function(measuArch,main){
#This R code allows us to reproduce the skeleton plots visualizing the seven archetypes of Figure 5 of the 
#paper Epifanio et al. (2013): "Archetypal analysis: Contributions for estimating boundary cases in 
#multivariate accommodation problem".

 #List with the measurements of each archetype:
 a <- measuArch
 
 #Popliteal height sitting:
 x1 <- c(40,40)
 y1 <- c(0,a[3]);
 par(pty = "s") 
 plot(x1,y1,type="l",lwd=5, xlim = c(0,60), ylim = c(0,60), xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      main = main)
 axis(1, at = seq(0,60,10), labels = seq(0,60,10))
 axis(2, at = seq(0,60,10), labels = seq(0,60,10))

 #Buttock knee length:
 x2 <- c(40, 40-a[2])
 y2 <- c(a[3], a[3])
 lines(x2,y2,lwd=5)

 #Shoulder height sitting:
 x3 <- c(40-a[2], 40-a[2])
 y3 <- c(a[3], a[3]+ a[6])
 lines(x3,y3,lwd=5)

 #Eye height sitting:
 x4 <- c(40-a[2], 40-a[2])
 y4 <- c(a[3], a[3]+ a[5])
 lines(x4,y4,lwd=5)

 #Sitting height:
 x5 <- c(40-a[2], 40-a[2])
 y5 <- c(a[3], a[3]+ a[4])
 lines(x5,y5,lwd=5)

 #Thumb tip reach:
 x6=c(40-a[2], 40-a[2]+a[1])
 y6=c(a[3]+a[6], a[3]+a[6])
 lines(x6,y6,lwd=5)

 #For the eyes postion:
 x7 <- c(40-a[2]-2, 40-a[2]+2)
 y7 <- c(a[3]+a[5], a[3]+a[5])
 lines(x7,y7,lwd=5)
}

