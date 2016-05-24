llbt.fit<-function(y,Xmodel,q,ncat,maxiter=100)
{
       Xmodel<-as.matrix(Xmodel)
       p<-ncol(Xmodel)
       if (qr(Xmodel)$rank<p)
          stop("linear dependencies in model matrix")
       b <- rep(0, p + q)
       dev <- 0
       fv <- y
       fv[y == 0] <- 1e-10
       eta <- log(fv)
       conv.dev.eps <- 1e-04
       iter <- 0
       #%ind.q <- 1:q
       #%Indmat.q <- cbind(ind.q, ind.q)
       repeat {
             b.old <- b
             dev.old <- dev
             z <- eta + (y - fv)/fv
             a22inv <- 1/colSums(matrix(fv, nrow = ncat))
             WX <- fv * cbind(z, Xmodel)
             A21 <- matrix(colSums(matrix(WX[, 1:(p + 1)], nrow = ncat)), ncol = p + 1)
             A11 <- crossprod(cbind(z, Xmodel), WX)
             A22inv.A21 <- -a22inv * A21
             A.11 <- solve(A11 + crossprod(A21, A22inv.A21))
             A.21 <- A22inv.A21 %*% A.11
             rm(A11, A21, a22inv, A22inv.A21, WX, z)
             pseudo.rss <- -1/A.11[1, 1]
             b.model <- pseudo.rss * A.11[1, 2:(p + 1)]
             b.elim <- pseudo.rss * A.21[, 1]
             b <- c(b.model, b.elim)
             #se <- sqrt(diag(A.11[2:(p + 1), 2:(p + 1)]) + pseudo.rss * A.11[2:(p + 1), 1] * A.11[2:(p + 1), 1]) #changed for case p==1
             se <- sqrt(A.11[cbind(2:(p + 1), 2:(p + 1))] + pseudo.rss * A.11[2:(p + 1), 1] * A.11[2:(p + 1), 1])
             rm(A.11, A.21)
             eta.model <- crossprod(t(Xmodel), b.model)
             eta.elim <- rep(b.elim, rep(ncat, q))
             eta <- as.vector(eta.model + eta.elim)
             fv <- exp(eta)
             dev <- sum(2 * (y * log(ifelse(y == 0, 1, y/fv)) - (y - fv)))
             iter <- iter + 1
             GLIM.dev.diff <- (dev - dev.old)/dev
             GLIM.converged.dev <- abs(GLIM.dev.diff) < conv.dev.eps
             if (iter > maxiter || GLIM.converged.dev) break
       }
       cat("\n")
       cat("Results of llbt.fit: \n")
       cat("\n")

       cat("Deviance:",dev,"\n")
       cat("Residual df =",length(y) - p - q,"\n")
       cat("Number of iterations:",iter,"\n")
       cat("\n")
       bb<-b[1:p]
       coeftable <- as.data.frame(cbind(round(bb,4),
                                  round(se,4),round(bb/se,4),round(1-pnorm(abs(bb/se)),4)))
       colnames(coeftable) <- c("Estimate","Std. Error","z","P(z)")
       if (is.null(colnames(Xmodel))){
          rownames(coeftable) <- paste("par",1:p,sep="")
       } else {
          rownames(coeftable) <-colnames(Xmodel)
       }
       print(coeftable)
}
