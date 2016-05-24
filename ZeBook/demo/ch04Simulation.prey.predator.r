# Francois Brun, 2013-03-05
# Predator-Prey Lotka-Volterra model (with logistic prey)
#### main program
library(deSolve)
library(ZeBook)
# method "euler" : Euler's Method,
# method "rk2" : Classical Runge-Kutta 2th Order Integration.
# method "rk4" : Classical Runge-Kutta 4th Order Integration.

system.time(sim <- predator.prey.model(grH=1, kH=10, mrH=0.2, eff=0.5, mrA=0.2, H0=1,A0=2,duration=200,dt=1, method="euler"))
plot(sim$time, sim$H, type="l", xlab = "time (day)", ylab = "population density",lty=1)
lines(sim$time, sim$A, lty=2)
legend("topright", c("H, prey", "A, predator"), lty = 1:2) 


# computation time for different integration methods and time step
seq_dt = c(0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
duration_euler = sapply(seq_dt, function(x) system.time(predator.prey.model(dt=x, method="euler"))["elapsed"]) 
duration_rk2 = sapply(seq_dt, function(x) system.time(predator.prey.model(dt=x, method="rk2"))["elapsed"]) 
duration_rk4 = sapply(seq_dt, function(x) system.time(predator.prey.model(dt=x, method="rk4"))["elapsed"]) 

plot(range(seq_dt),c(0,max(duration_rk4)), xlab="dt", ylab="computer time (s)", type="n", log="x")
lines(seq_dt,duration_euler, lty=1, lwd=2)
lines(seq_dt,duration_rk2, lty=2, lwd=2)
lines(seq_dt,duration_rk4, lty=3, lwd=2)
legend("topright", c("Euler", "rk2", "rk4"), lty = 1:3,lwd=2) 

# end of file