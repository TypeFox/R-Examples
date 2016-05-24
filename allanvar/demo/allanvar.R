#Loading values
data(gyroz)

#Allan variance computation using avar
avgyroz <- avar(gyroz@.Data[1:1000], frequency(gyroz))
plotav(avgyroz)
abline(1.0+log(avgyroz$time[1]), -1/2, col="green", lwd=4, lty=10)
abline(1.0+log(avgyroz$time[1]), 1/2, col="green", lwd=4, lty=10)
legend(0.11, 1e-03, c("Random Walk"))
legend(2, 1e-03, c("Rate Random Walk"))

#Allan variance computation using avarn
avngyroz <- avarn(gyroz@.Data[1:1000], frequency(gyroz))
plotav(avngyroz)
abline(1.0+log(avngyroz$time[1]), -1/2, col="green", lwd=4, lty=10)
abline(1.0+log(avngyroz$time[1]), 1/2, col="green", lwd=4, lty=10)
legend(0.11, 1e-03, c("Random Walk"))
legend(2, 1e-03, c("Rate Random Walk"))

##Allan variance computation using avari
##Simple integration of the angular velocity
igyroz <- cumsum(gyroz@.Data[1:1000] * (1/frequency(gyroz)))
igyroz <- ts (igyroz, start=c(igyroz[1]), delta=(1/frequency(gyroz)))
avigyroz <- avari(igyroz@.Data, frequency(igyroz))
plotav(avigyroz)
abline(1.0+log(avigyroz$time[1]), -1/2, col="green", lwd=4, lty=10)
abline(1.0+log(avigyroz$time[1]), 1/2, col="green", lwd=4, lty=10)
legend(0.11, 1e-03, c("Random Walk"))
legend(2, 1e-03, c("Rate Random Walk"))


#Ploting all
plot (avgyroz$time,sqrt(avgyroz$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avngyroz$time,sqrt(avngyroz$av), col="green", lwd=1)
lines (avigyroz$time,sqrt(avigyroz$av), col="red")
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis Comparison", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroZ (avar)", "GyroZ(avarn)", "GyroZ(avari)"),  fill = c("blue", "green", "red"))
