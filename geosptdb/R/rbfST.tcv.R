assign("rbfST.tcv",
       function (formula, data, eta, rho, n.neigh, func, progress=TRUE) 
       {
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
         if(progress)
         pb <- txtProgressBar(min = 0, max = length(z), char = "=", style = 3)
         for(i in 1:(length(z))){                                                
           pred[i] <-  rbf0(S=S[-i,], z, X0=as.matrix(X0[-i,]), x0=as.matrix(X0[i,]), dist.S=dist.S[i,-i], eta, rho, n.neigh, func)  #  s=S[i,]
           if(progress)
           setTxtProgressBar(pb, i)
         }
         if(progress)
         close(pb)
         rbf.pred <- data.frame(pred,NA,z,NA,NA,NA,s0)
         colnames(rbf.pred) <- c("var1.pred","var1.var","observed","residual","zscore","fold","x","y","t")
         rbf.pred[,6] <- 1:length(z)
         rbf.pred[,3]<- z
         rbf.pred[,4]<- rbf.pred[,3]-rbf.pred[,1]
         rbf.pred
       }
)
