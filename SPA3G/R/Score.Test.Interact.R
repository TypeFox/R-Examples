Score.Test.Interact <-
function(Y, K1, K2, K3, par, method="BFGS", test=TRUE)
{
n <- length(Y)
p <- length(par)
if (p!=3 & test==TRUE)
cat("Error: Not matched initial values!")
theta.new <- par
theta.old <- rep(0, p)
     X <- matrix(1, n, 1)
Vs <- array(0, c(n, n, 4))
Vs[, , 1] <- diag(1, n)
Vs[, , 2] <- K1
Vs[, , 3] <- K2
Vs[, , 4] <- K3
Sigma <- 0
for (i in 1 : p)
{
Sigma <- Sigma+theta.new[i]*Vs[, , i]
}
W <- solve(Sigma)  
R <- W-W%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W
kk <- g.old <- 0; tt <- c()
#lR <- 1; lR.old <- 0
while(sum(abs(theta.new-theta.old))>1.0e-05 & kk <100)
#while(abs(lR-lR.old)>1.0e-05 & kk <100)
{
#lR.old <- lR
if (method=="BFGS")
{
s <- theta.new-theta.old
theta.old <- theta.new
g <- c()
for (i in 1 : p)
{
              g[i] <- -t(Y)%*%R%*%Vs[, , i]%*%R%*%Y+TT(R, Vs[, , i])
            }
delta <- g-g.old
g.old <- g
if (kk==0 | t(s)%*%delta <=0)
{
AI <- matrix(0, p, p)
for (i in 1 : p)
{
for (j in i : p)
{
AI[i, j] <- AI[j, i] <- t(Y)%*%R%*%Vs[, , i]%*%R%*%Vs[, , j]%*%R%*%Y
                }
            }
H_inv <- solve(AI)
} else
{
rho <- c(1/(t(delta)%*%s))
H_inv <- (diag(1, p)-(s%*%t(delta))*rho)%*%H_inv%*%(diag(1, p)-rho*delta%*%t(s))+rho*s%*%t(s)
}
}
if (method=="AI")
{
theta.old <- theta.new
g <- c()
for (i in 1 : p)
{
              g[i] <- t(Y)%*%R%*%Vs[, , i]%*%R%*%Y-TT(R, Vs[, , i])
            }
H <- matrix(0, p, p)
for (i in 1 : p)
{
for (j in i : p)
{
H[i, j] <- H[j, i] <- -t(Y)%*%R%*%Vs[, , i]%*%R%*%Vs[, , j]%*%R%*%Y
                }
            }
H_inv <- solve(H)
}

if (method=="FS")
{
theta.old <- theta.new
g <- c()
for (i in 1 : p)
{
              g[i] <- t(Y)%*%R%*%Vs[, , i]%*%R%*%Y-TT(R, Vs[, , i])
          }
H <- matrix(0, p, p)
for (i in 1 : p)
{
AA <- R%*%Vs[, , i]
for (j in i : p)
{
BB <- R%*%Vs[, , j]
H[i, j] <- H[j, i] <- -TRACE(AA%*%BB) 
                }
            }
H_inv <- solve(H)
}

theta.new <- theta.old-H_inv%*%(g)
alpha <- 0.5
while (length(which(theta.new<0))>0 & alpha>1.0e-008)
{
theta.new <- theta.old-alpha*H_inv%*%(g)
            alpha <- alpha/2
      }
theta.new[which(theta.new<0)] <- 0
Sigma.new <- 0
for (i in 1 : p)
{
Sigma.new <- Sigma.new+theta.new[i]*Vs[, , i]
}
W.new <- solve(Sigma.new)
R <- W.new-W.new%*%X%*%solve(t(X)%*%W.new%*%X)%*%t(X)%*%W.new
#lR <- -(log(det(Sigma.new))+log(det(t(X)%*%W.new%*%X))+t(Y)%*%R%*%Y)/2
#tt <- c(tt, lR)
#plot(tt)
kk <- kk+1
}

a1 <- R%*%Vs[, , 1]
a2 <- R%*%Vs[, , 2]
a3 <- R%*%Vs[, , 3]
a4 <- R%*%Vs[, , 4]
b11 <- TT(a1, a1)
b12 <- TT(a1, a2)
b13 <- TT(a1, a3)
b14 <- TT(a1, a4)
b22 <- TT(a2, a2)
b23 <- TT(a2, a3)
b24 <- TT(a2, a4)
b33 <- TT(a3, a3)
b34 <- TT(a3, a4)
b44 <- TT(a4, a4)
 
if (test==FALSE)
{
eigen.sigma <- eigen(Sigma.new)
lR <- -(sum(log(eigen.sigma$values))+log(det(t(X)%*%W.new%*%X))+t(Y)%*%R%*%Y)/2
H <- matrix(c(b11, b12, b13, b14, b12, b22, b23, b24, b13, b23, b33, b34, b14, b24, b34, b44),4, 4)/2
beta <- solve(t(X)%*%W.new%*%X)%*%t(X)%*%W.new%*%Y
object <- list(VCs=theta.new, fisher.info=H, Beta=beta, restricted.logLik=lR)
return(object)
}

if (test==TRUE)
{
eigen.sigma <- eigen(Sigma.new)
lR <- -(sum(log(eigen.sigma$values))+log(det(t(X)%*%W.new%*%X))+t(Y)%*%R%*%Y)/2
W0 <- W.new
beta <- solve(t(X)%*%W0%*%X)%*%t(X)%*%W0%*%Y
Q <- t(Y-X%*%beta)%*%W0%*%K3%*%W0%*%(Y-X%*%beta)/2
e <- TT(R, K3)/2
Its <- c(b14, b24, b34)
Iss <- matrix(c(b11, b12, b13, b12, b22, b23, b13, b23, b33),3, 3)
Itt <- (b44-Its%*%solve(Iss)%*%Its)/2
k <- Itt/e/2; v=2*e^2/Itt
pvalue <- pchisq(Q/k, df=v, lower.tail=F)
object <- list(VCs=theta.new, fisher.info=Iss/2, Beta=beta, restricted.logLik=lR, 
   Score=Q, df=v, scale=k, p.value=pvalue)
class(object) <- "Score Test: tau3=0"
return(object)
} 
}
