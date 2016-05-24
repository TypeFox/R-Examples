assign("rbf",
       function (formula, data, eta, rho, newdata, n.neigh, func)
       {
         s = coordinates(data)
         s0 = coordinates(newdata)
         z = extractFormula(formula, data, newdata)$z
         X0 = extractFormula(formula, data, newdata)$X
         x0 = extractFormula(formula, data, newdata)$x0
         dist.newdata <- rdist(s, s0)
         vec.orden <- apply(dist.newdata, 2, function(x) x <- x[order(x)][1:n.neigh])
         vec.orden1 <- apply(dist.newdata, 2, function(x) x <- order(x)[1:n.neigh])
         Dist <- lapply(1:nrow(coordinates(newdata)), function(i) rdist(s[vec.orden1[,i],]))
         I.PHI.M <- lapply(1:nrow(coordinates(newdata)), function(i) {
           if (func == "M" | func == "ST" | func == "CRS" | func == "TPS" | func == "TRI")
             solve(RBF.phi(Dist[[i]] - diag(Dist[[i]]), eta, func) + rho * diag(n.neigh))
           else chol2inv(chol(RBF.phi(Dist[[i]] - diag(Dist[[i]]), eta, func) + rho * diag(n.neigh)))
         })
         b <- RBF.phi(vec.orden, eta, func)
         X <- vector(mode = "list", length = nrow(coordinates(newdata)))
         Xi <- lapply(X, function(x) x <- X0)
         XT <- lapply(1:nrow(coordinates(newdata)), function(i) Xi[[i]][vec.orden1[,i],])
         rm(vec.orden, Xi)
         Lambda <- lapply(1:nrow(coordinates(newdata)), function(i) 
         I.PHI.M[[i]] %*%(b[, i] - as.matrix(XT[[i]]) %*% (Solve(t(XT[[i]]) %*% I.PHI.M[[i]] %*% XT[[i]])) %*% 
                              (t(XT[[i]]) %*% I.PHI.M[[i]] %*% b[,i] - x0[i, ])))
         remove(I.PHI.M, XT, x0)
         pred <- lapply(1:nrow(coordinates(newdata)), function(i) t(Lambda[[i]])%*%z[vec.orden1[, i]])
         rm(vec.orden1)
         rbf.pred <- as.data.frame(matrix(NA, nrow = nrow(coordinates(newdata)),ncol = 4))
         colnames(rbf.pred) <- c("x", "y", "var1.pred", "var1.var")
         rbf.pred[, 1:2] <- s0
         rbf.pred[, 3] <- unlist(pred)
         rbf.pred
       }
)
