assign("rbfST.cv1",  
       function (param, formula, data, n.neigh, func) 
       {
eta <- param[1]
rho <- param[2]
s = cbind(coordinates(data),data["t"]@data)
s0 = cbind(coordinates(data),data["t"]@data)
z = extractFormula(formula, data, data)$z
X0 = extractFormula(formula, data, data)$X
S <- scale(s)  
dist.S <- rdist(S)-diag(rdist(S))
rbf0 <- function(S, z, X0, x0, dist.S, eta, rho, n.neigh, func){     
  vec.orden <- order(dist.S)
  vc <- vec.orden[1:n.neigh]
  dist.vec.cerca <- dist.S[vc]
  X <- X0[vc,]
  M.dist <- rdist(S[vc,])
  m.dist.vec <- M.dist-diag(M.dist)
  phi <- RBF.phi(m.dist.vec,eta,func)
  PHI <- phi+rho*diag(n.neigh)
  b <- RBF.phi(dist.vec.cerca,eta,func)
  I.PHI <- if (func %in% c("M","ST","CRS","TPS","TRI")) solve(PHI)
  else chol2inv(chol(PHI))
  Lambda <- I.PHI%*%(b-as.matrix(X)%*%(Solve(t(X)%*%I.PHI%*%X))%*%(t(X)%*%I.PHI%*%b-x0))
  pred <- t(Lambda)%*%z[vc]
  pred
}      
pred <- as.numeric(NA,length= length(z))
#pb <- txtProgressBar(min = 0, max = length(z), char = "=", style = 3)
for(i in 1:(length(z))){                                                
  pred[i] <-  rbf0(S=S[-i,], z, X0=as.matrix(X0[-i,]), x0=as.matrix(X0[i,]), dist.S=dist.S[i,-i], eta, rho, n.neigh, func)  
#  setTxtProgressBar(pb, i)
}
#close(pb)
rbf.pred <- data.frame(s0,pred,NA)
colnames(rbf.pred) <- c("x","y","t","var1.pred","var1.var")
RMSPE <-  sqrt(sum((rbf.pred$var1.pred-z)^2)/length(z))
RMSPE
}
)

