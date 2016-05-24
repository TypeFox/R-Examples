###################################################
### chunk number 1: 
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### chunk number 2: 
###################################################
exp2d.data <- exp2d.rand(n2=150, lh=0, dopt=10)
X <- exp2d.data$X
Z <- exp2d.data$Z
XX <- rbind(c(0,0),c(2,2),c(4,4))


###################################################
### chunk number 3: 
###################################################
out <- btgpllm(X=X, Z=Z, XX=XX, corr="exp", bprior="b0", 
               pred.n=FALSE, Ds2x=TRUE, R=10, #BTE=c(2000,5000,10),
               trace=TRUE, verb=0)


###################################################
### chunk number 4: 
###################################################
out$trace


###################################################
### chunk number 5: XXd
###################################################
trXX <- out$trace$XX; ltrXX <- length(trXX)
y <- trXX[[1]]$d
for(i in 2:ltrXX) y <- c(y, trXX[[i]]$d)
plot(log(trXX[[1]]$d), type="l", ylim=range(log(y)), ylab="log(d)",
     main="range (d) parameter traces")
names <- "XX[1,]"
for(i in 2:ltrXX) {
  lines(log(trXX[[i]]$d), col=i, lty=i)
  names <- c(names, paste("XX[", i, ",]", sep=""))
}
legend("bottomleft", names, col=1:ltrXX, lty=1:ltrXX)


###################################################
### chunk number 6: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 7: 
###################################################
linarea <- mean(out$trace$linarea$la)
linarea


###################################################
### chunk number 8: la
###################################################
hist(out$trace$linarea$la)


###################################################
### chunk number 9: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 10: 
###################################################
m <- matrix(0, nrow=length(trXX), ncol=3)#ncol=5)
for(i in 1:length(trXX))
  m[i,] <- as.double(c(out$XX[i,], mean(trXX[[i]]$b)))
m <- data.frame(cbind(m, 1-m[,3]))
names(m)=c("XX1","XX2","b","pllm")
m


###################################################
### chunk number 11: alc
###################################################
trALC <- out$trace$preds$Ds2x
y <- trALC[,1]
for(i in 2:ncol(trALC)) y <- c(y, trALC[,i])
plot(log(trALC[,1]), type="l", ylim=range(log(y)), ylab="Ds2x",
     main="ALC: samples from Ds2x")
names <- "XX[1,]"
for(i in 2:ncol(trALC)) {
  lines(log(trALC[,i]), col=i, lty=i)
  names <- c(names, paste("XX[", i, ",]", sep=""))
}
legend("bottomright", names, col=1:ltrXX, lty=1:ltrXX)


###################################################
### chunk number 12: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


