slcestimates <- function(info,n_a)  {
slength <- length(info)
n_b <- slength-n_a
phaseA <- info[1:n_a]
phaseB <- info[(n_a+1):slength]
# Estimate phase A trend
phaseAdiff <- c(1:(n_a-1))
for (iter in 1:(n_a-1))
phaseAdiff[iter] <- phaseA[iter+1] - phaseA[iter]
trendA <- mean(phaseAdiff)
# Remove phase A trend from the whole data series
phaseAdet <- c(1:n_a)
for (timeA in 1:n_a)
phaseAdet[timeA] <- phaseA[timeA] - trendA * timeA
phaseBdet <- c(1:n_b)
for (timeB in 1:n_b)
phaseBdet[timeB] <- phaseB[timeB] - trendA * (timeB+n_a)
# Compute the slope change estimate
phaseBdiff <- c(1:(n_b-1))
for (iter in 1:(n_b-1))
phaseBdiff[iter] <- phaseBdet[iter+1] - phaseBdet[iter]
trendB <- mean(phaseBdiff)
print ("Slope change estimate = "); print(trendB)
# Compute the level change estimate
phaseBddet <- c(1:n_b)
for (timeB in 1:n_b)
phaseBddet[timeB] <- phaseBdet[timeB] - trendB * (timeB-1)
level <- mean(phaseBddet) - mean(phaseAdet)
print ("Level change estimate = "); print(level)

# Represent graphically
time <- c(1:slength)
par(mfrow=c(2,1))
plot(time,info, xlim=c(1,slength), ylim=c((min(info)-1),(max(info)+1)), xlab="Measurement time", ylab="Variable of interest", font.lab=2)
abline(v=(n_a+0.5))
lines(time[1:n_a],info[1:n_a])
lines(time[(n_a+1):slength],info[(n_a+1):slength])
axis(side=1, at=seq(0,slength,1),labels=TRUE, font=2)
axis(side=2, at=seq((min(info)-1),(max(info)+1),2),labels=TRUE, font=2)
points(time, info, pch=24, bg="black")
title (main="Original data")

transf <- c(phaseAdet,phaseBdet)
plot(time,transf, xlim=c(1,slength), ylim=c((min(transf)-1),(max(transf)+1)), xlab="Measurement time", ylab="Variable of interest", font.lab=2)
abline(v=(n_a+0.5))
lines(time[1:n_a],transf[1:n_a])
lines(time[(n_a+1):slength],transf[(n_a+1):slength])
axis(side=1, at=seq(0,slength,1),labels=TRUE, font=2)
axis(side=2, at=seq((min(transf)-1),(max(transf)+1),2),labels=TRUE, font=2)
points(time, transf, pch=24, bg="black")
title (main="Detrended data")

list(trendB,level)

}
