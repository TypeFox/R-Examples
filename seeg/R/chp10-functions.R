
lda2 <- function(X1,X2){

# two-group linear discriminant analysis
# Miguel Acevedo march 2004 

# number of variables and observations in each group
m <-  dim(X1)[2]
n1 <- dim(X1)[1]
n2 <- dim(X2)[1]

# calculate sample means 
G1 <- seq(1,m); G2 <- G1
for(i in 1:m) G1[i] <- mean(X1[,i])
for(i in 1:m) G2[i] <- mean(X2[,i])
# difference of sample means
D <- G1 - G2
# var/cov
S1 <- var(X1)
S2 <- var(X2)
# pooled var-cov
Sp <- (1/(n1+n2-2)) * ((n1-1)*S1 + (n2-1)*S2)
# invert
I <- diag(m)
Sinv <- solve(Sp,I)

# coefficients
A <- Sinv%*%D
# scores of means
Z1 <- G1%*%A 
Z2 <- G2%*%A
# score midpoint 
Z0 <- (1/2)*(G1+G2)%*%A

# Mahalanobis distance
D2 <- t(D)%*%Sinv%*%D
# equivalently
#D2 <- t(D)%*%A
#D2 <- abs(Z1-Z2)

# observation scores centered
Z1s <- c(X1%*%A) - Z0
Z2s <- c(X2%*%A) - Z0

# centroid scores centered 
Z1c <- Z1 - Z0
Z2c <- Z2 - Z0

smin <- min(Z1s,Z2s); smax <- max(Z1s,Z2s)
# plot
panel3(size=7)
plot(Z1s,rnorm(n1), xlim=c(smin,smax),ylim=c(-10,10), col=1, 
     ylab="Arbitrary",xlab="Discriminant Scores")
points(Z2s,rnorm(n2), pch=2, col=1)
abline(h=0)
points(Z1c,0, pch="1", col=1,cex=2)
points(Z2c,0, pch="2", col=1,cex=2)

hist(Z1s,xlim=c(smin,smax),prob=T)
hist(Z2s,xlim=c(smin,smax),prob=T)

#Two-sample version of Hotelling’s T2 test (similar to two sample t test)
 dof.num <- m
 dof.den <- n1 + n2 - m + 1
 if(dof.den > 0) {
  F1 <- (n1 +n2 -m -1)/ (m*(n1+n2-2))
  F2 <- (n1*n2)/(n1+n2)
  F <- F1*F2*D2
  p.value <- 1 - pf(F, dof.num, dof.den)
  }
 else {
  D2 <- NA; F <- NA; p.value <- NA
  }
m.n1.n2 <- c(m,n1,n2)
Z1c.Z2c <-  c(Z1c, Z2c)
D2.F.p.value <- c(D2, F, p.value)

return(list(m.n1.n2=m.n1.n2, G1=G1, G2=G2, S1=S1, S2=S2, Sp=Sp, Spinv=Sinv, A=A, 
       Z1c.Z2c=Z1c.Z2c, D2.F.p.value=D2.F.p.value, Z1s=Z1s, Z2s=Z2s))

} 	
