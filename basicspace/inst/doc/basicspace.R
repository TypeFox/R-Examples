### R code from vignette source 'basicspace.Rnw'

###################################################
### code chunk number 1: one
###################################################
set.seed(1231)
library("basicspace")
N <- 1000
M <- 20
s <- 2
fraction.missing <- 0.3
E <- matrix(runif(N * M, min = -0.5, max = 0.5), nrow = N, ncol = M)


###################################################
### code chunk number 2: two
###################################################
U <- matrix(runif(N * s), nrow = N, ncol = s)
D <- diag(seq(from = 2.1, by = -0.2, length.out = s))
V.prime <- matrix(runif(s * M), nrow = s, ncol = M)
c <- rnorm(M)
Jn <- rep(1, N)


###################################################
### code chunk number 3: three
###################################################
X.true <- U %*% D %*% V.prime + Jn %o% c
X.0 <- X.true + E
Psi.true <- U %*% sqrt(D)
W.true <- t(V.prime) %*% sqrt(D)


###################################################
### code chunk number 4: four
###################################################
missing <- sample(1:(N * M), round(fraction.missing * N * M))
X.0[missing] <- 999


###################################################
### code chunk number 5: five
###################################################
rownames(X.0) <- paste("Legis", 1:N, sep = "")
colnames(X.0) <- paste("V", 1:M, sep = "")


###################################################
### code chunk number 6: six
###################################################
result <- blackbox(X.0, missing = c(999), verbose = TRUE, dims = 3, minscale = 8)
names(result)


###################################################
### code chunk number 7: seven
###################################################

Psi.hat <- cbind(result$individuals[[2]]$c1, result$individuals[[2]]$c2)
c.hat <- result$stimuli[[2]]$c
xrow <- sapply(1:N, function(i) length(rep(1,s)[!is.na(Psi.hat[i,])]))
Psi.hat <- Psi.hat[!(xrow<2),]
Psi.true <- Psi.true[!(xrow<2),]
Psi.hat[,1] <- Psi.hat[,1]-mean(Psi.hat[,1])
Psi.hat[,2] <- Psi.hat[,2]-mean(Psi.hat[,2])
Psi.true[,1] <- Psi.true[,1]-mean(Psi.true[,1])
Psi.true[,2] <- Psi.true[,2]-mean(Psi.true[,2])
C <- t(Psi.true)%*%Psi.hat
svddecomp <- svd(C)
U.rotate <- svddecomp$u
V.rotate <- svddecomp$v
T <- V.rotate %*% t(U.rotate)
Psi.hatrotate <- Psi.hat %*% T

par(mfrow=c(1,2))
plot(Psi.true[,1],Psi.hatrotate[,1], xlim= c(-0.7,0.7), ylim= c(-0.5, 0.5),
pch=20, cex=0.4, cex.lab=1.6, bty="n",
xlab="True Psi, first dimension",
ylab="Recovered Psi, first dimension")
plot(Psi.true[,2],Psi.hatrotate[,2], xlim= c(-0.7,0.7), ylim= c(-0.5, 0.5),
pch=20, cex=0.4, cex.lab=1.6, bty="n",
xlab="True Psi, second dimension", ylab="Recovered Psi, second dimension")


###################################################
### code chunk number 8: six
###################################################
W.hat <- cbind(result$stimuli[[2]]$w1, result$stimuli[[2]]$w2)
W.hatrotate <- W.hat%*%T
par(mfrow=c(1,2))

plot(W.true[,1],W.hatrotate[,1], xlim= c( 0.00, 1.50), ylim= c(-0.75, 2.75),
pch=20, cex=1.5, cex.lab=1.6, bty="n",
xlab="True W, first dimension", ylab="Recovered W, first dimension")

plot(W.true[,2],W.hatrotate[,2], xlim= c( 0.00, 1.50), ylim= c(-0.75, 2.75),
pch=20, cex=1.5, cex.lab=1.6, bty="n",
xlab="True W, second dimension", ylab="Recovered W, second dimension")



###################################################
### code chunk number 9: eight
###################################################

par(mfrow=c(1,1))

plot(c, c.hat, 
pch=20, cex=1.2, cex.lab=1.1, bty="n", 
xlab="True C", ylab="Recovered C")


###################################################
### code chunk number 10: nine
###################################################

W.hat <- cbind(result$stimuli[[2]]$w1, result$stimuli[[2]]$w2)
Psi.hat <- cbind(result$individuals[[2]]$c1, result$individuals[[2]]$c2)
X.hat <- Psi.hat %*% t(W.hat) + Jn %o% result$stimuli[[2]]$c

par(mfrow=c(1,2))

plot(X.true[missing], X.hat[missing],
pch = 20, cex = 0.4, cex.lab = 1.2, bty = "n",
xlab = "True X, missing values", ylab = "Recovered X, missing values")

plot(X.true[!(1:(N*M) %in% missing)], X.hat[!(1:(N*M) %in% missing)],
pch = 20, cex = 0.4, cex.lab = 1.2, bty = "n",
xlab = "True X, nonmissing values", ylab = "Recovered X, nonmissing values")


###################################################
### code chunk number 11: ten
###################################################
data("Issues1980")
Issues1980[1:10, 1:4]


###################################################
### code chunk number 12: eleven
###################################################
Issues1980[Issues1980[, "abortion1"] == 7, "abortion1"] <- 8
Issues1980[Issues1980[, "abortion2" ]== 7, "abortion2"] <- 8


###################################################
### code chunk number 13: twelve
###################################################
result <- blackbox(Issues1980, missing = c(0,8,9), verbose = FALSE, dims = 2,
   minscale = 8)


###################################################
### code chunk number 14: thirteen
###################################################
summary(result)


###################################################
### code chunk number 15: fourteen
###################################################
cor(result$individuals[[1]]$c1, Issues1980[, "libcon1"], use="pairwise")


###################################################
### code chunk number 16: blackbt
###################################################
data("LC1980")
LCdat=LC1980[, -1]
LCdat[1:10,]
result <- blackbox_transpose(LCdat, missing = c(0,8,9), dims = 2,
   minscale = 5, verbose = TRUE)


###################################################
### code chunk number 17: bbt
###################################################

par(mfrow = c(1, 2))

plot(result)
plotcdf.blackbt(result)


###################################################
### code chunk number 18: blackfinal
###################################################
summary(result)


###################################################
### code chunk number 19: aldmck
###################################################
data("LC1980")
result <- aldmck(data = LC1980, polarity = 2, respondent = 1,
    missing = c(0, 8, 9), verbose = TRUE)
summary(result)


###################################################
### code chunk number 20: aldmck_plot
###################################################

plot.aldmck(result)


###################################################
### code chunk number 21: aldmck_boot
###################################################
result <- boot_aldmck(data=LC1980, polarity=2, respondent=1,
    missing=c(0,8,9), iter=100)
apply(result, 2, sd)


###################################################
### code chunk number 22: aldmck_MC
###################################################
Nstimuli <- 6
Nresp <- 500

Z_j <- rnorm(6)
Z_j <- (Z_j - mean(Z_j))/sd(Z_j)

respondent.sd <- runif(Nresp, min = 0.3, max = 0.9)
error_heteroskedastic <- matrix(NA, Nresp, Nstimuli)
for(i in 1:Nresp) error_heteroskedastic <- rnorm(Nstimuli, sd = respondent.sd)

w_i <- runif(Nresp, min=0, max=1)
c_i <- rnorm(Nresp)
Y_ij <- rep(1,500) %o%  Z_j
Y_ij <- Y_ij + error_heteroskedastic
R_ij <-  1/w_i %o% rep(1,Nstimuli) * (Y_ij - c_i %o% rep(1,Nstimuli))

result <- aldmck(R_ij, polarity = 6, missing = c(999))

cor(Z_j, result$stimuli)


