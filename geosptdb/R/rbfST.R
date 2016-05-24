assign("rbfST",
       function(formula, data, eta, rho, newdata, n.neigh, func, progress=TRUE){
         s = cbind(coordinates(data),data["t"]@data)
         s0 = cbind(coordinates(newdata),newdata["t"]@data)
         z = extractFormula(formula, data, newdata)$z
         X0 = extractFormula(formula, data, newdata)$X
         x0 = extractFormula(formula, data, newdata)$x0
         S <- scale(s)
         dist.newdata <- rdist(S, cbind(standardize(s0[,1],mean(s[,1]),sd(s[,1])),standardize(s0[,2],mean(s[,2]),sd(s[,2])),
                                        standardize(s0[,3],mean(s[,3]),sd(s[,3]))))
         Pred <- as.numeric(NA,length= nrow(coordinates(newdata)))
         rbf0 <- function(S, z, X0, x0, dist.newdata, eta, rho, n.neigh, func){
           vec.orden <- order(dist.newdata)
           vc <- vec.orden[1:n.neigh]
           dist.vec.cerca <- dist.newdata[vc]
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
         if(progress)
         pb <- txtProgressBar(min = 0, max = nrow(coordinates(newdata)), char = "=", style = 3)
         for (i in 1:nrow(coordinates(newdata))){
           Pred[i] <- rbf0(S, z, X0, x0=x0[i,], dist.newdata=dist.newdata[,i], eta, rho, n.neigh, func)
           if(progress)
           setTxtProgressBar(pb, i)
         }
         if(progress) 
         close(pb)
         rbf.pred <- data.frame(s0,Pred,NA)
         names(rbf.pred) <- c("x","y","t","var1.pred","var1.var")
         rbf.pred
       }
)  
