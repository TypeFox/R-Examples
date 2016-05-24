emul.predict <-
function(emul, theta.star) {

#AN EXCEPTION FOR EXTRAPOLATION #!+
par.min             <- apply(emul$Theta.mat, 2, min)
par.max             <- apply(emul$Theta.mat, 2, max)
anyhigh             <- any(theta.star > par.max)
anylow              <- any(theta.star < par.min)
if (anyhigh || anylow ) stop("***ERROR*** Prediction point is outside of ensemble range")


#PRELIMINARIES #!+
I.mat               <- diag(emul$n)#!+
Sigma.mats          <- sep.cov(emul$Theta.mat, emul$t.vec, emul$rho, emul$kappa,
                      emul$phi.vec, emul$zeta) #!+
Sigma.theta.inv.mat <- solve(Sigma.mats$Sigma.theta.mat) #!+
n.par               <- emul$n #!+
p.par               <- emul$p #!+
m.par               <- dim(emul$Theta.mat)[2] #!+
Theta.mat           <- emul$Theta.mat #!+

# MU VECTOR #!+
# The first term in Equation 11
X.star.mat          <- matrix(1, nrow=n.par, ncol=1) 
Theta.star.mat      <- matrix(theta.star, nrow=n.par, ncol=m.par, byrow=TRUE)
if (any(emul$par.reg)) { 
  X.star.mat        <- cbind(X.star.mat, Theta.star.mat[,emul$par.reg])
}
if (emul$time.reg)     X.star.mat  <- cbind(X.star.mat, emul$t.vec) 
mu.vec              <- X.star.mat%*%as.matrix(emul$beta.vec)


# CONSTRUCT SIGMA.TH.STAR.TH MATRIX #!+
Theta.star.pmat     <- matrix(theta.star, nrow=p.par, ncol=m.par, byrow=TRUE)
Phi.mat             <- matrix(emul$phi.vec, nrow=p.par, ncol=m.par, byrow=TRUE)
Phisq.mat           <- Phi.mat^2 
Diff.mat            <- ((Theta.star.pmat - Theta.mat)^2)/Phisq.mat 
M1                  <- t(apply(Diff.mat, 1, sum)) 
Sigma.th.star.th    <- emul$kappa*exp(-M1) 
same.vec            <- M1 == 0 # #Nonzero elements indicate that Theta*=Theta_j
Sigma.th.star.th[same.vec] <- Sigma.th.star.th[same.vec] + emul$zeta 


# PREDICTIVE MEAN #w
M2             <- Sigma.th.star.th%*%Sigma.theta.inv.mat #w

#Simplified computation of M3 matrix without the kronecker product #!+
M3 <- matrix(0, nrow=n.par, ncol=n.par*p.par) #!+
# Fill M3 matrix by row #!+
for (myrow in 1:n.par) {
   start.col <- p.par*(myrow-1) +1 
   end.col   <- myrow*p.par        
   M3[myrow,start.col:end.col] <- M2 
}
T2             <- M3%*%emul$vecC #w
mu.star.vec    <- mu.vec + T2#w

# PREDICTIVE STANDARD DEVIATION #!+
# It is forced to zero if we are predicting at one of the design points
Term1          <- (emul$kappa + emul$zeta)*Sigma.mats$Sigma.t.mat 
M4             <- Sigma.th.star.th%*%Sigma.theta.inv.mat%*%t(Sigma.th.star.th)
Term2          <- as.vector(M4)*Sigma.mats$Sigma.t.mat #!+
Sigma.star.mat <- Term1 - Term2 
if (any(same.vec)) { 
  Sigma.star.mat <- matrix(0, nrow=n.par, ncol=n.par)
}

# FORMAT OUTPUT #!+
predict.out    <- list(mean=mu.star.vec, covariance=Sigma.star.mat) 
predict.out
}
