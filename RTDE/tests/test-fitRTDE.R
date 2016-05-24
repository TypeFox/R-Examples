
library(RTDE)


#####
# (1) simulation 

omega <- 1/2
m <- 48
n <- 100
obs <- cbind(rupareto(n), rupareto(n)) + rupareto(n)

#function of m
system.time(
x <- fitRTDE(obs, nbpoint=m:(n-m), 0, 1/2)
)
x
summary(x)
plot(x, which=1)
plot(x, which=2)

#function of m, alpha
system.time(
x <- fitRTDE(obs, nbpoint=m:(n-m), alpha=c(0, .25), omega=omega, control=list(trace=0))
)
x
dim(x$eta)
dimnames(x$eta)
summary(x)
str(x)

plot(x, which=1)
plot(x, which=2)


#function of m, alpha, omega
system.time(
x <- fitRTDE(obs, nbpoint=m:(n-m), alpha=c(0, .25), omega=c(1/2, 1/3), control=list(trace=0))
)
x
dim(x$eta)
dimnames(x$eta)
summary(x)

plot(x, which=1)
plot(x, which=2)
