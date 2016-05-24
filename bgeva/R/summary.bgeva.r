summary.bgeva <- function(object,s.meth="svd",sig.lev=0.05,...){

  testStat <- function (p, X, V, rank = NULL) {
      qrx <- qr(X)
      R <- qr.R(qrx)
      V <- R %*% tcrossprod(V[qrx$pivot, qrx$pivot],R)
      V <- (V + t(V))/2
      ed <- eigen(V, symmetric = TRUE)
      k <- max(0, floor(rank))
      nu <- abs(rank - k)

          if (rank > k + 0.05 || k == 0) 
              k <- k + 1
          nu <- 0
          rank <- k
      
      if (nu > 0) 
          k1 <- k + 1
      else k1 <- k
      r.est <- sum(ed$values > max(ed$values) * .Machine$double.eps^0.9)
      if (r.est < k1) {
          k1 <- k <- r.est
          nu <- 0
          rank <- r.est
      }
      vec <- ed$vectors
      if (k1 < ncol(vec)) 
          vec <- vec[, 1:k1, drop = FALSE]
      if (k == 0) {
          vec <- t(t(vec) * sqrt(nu/ed$val[1]))
      }
      if (nu > 0 && k > 0) {
          if (k > 1) 
              vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k - 
                  1)]))
          b12 <- 0.5 * nu * (1 - nu)
          if (b12 < 0) 
              b12 <- 0
          b12 <- sqrt(b12)
          B <- matrix(c(1, b12, b12, nu), 2, 2)
          ev <- diag(ed$values[k:k1]^-0.5)
          B <- ev %*% B %*% ev
          eb <- eigen(B, symmetric = TRUE)
          rB <- eb$vectors %*% tcrossprod(diag(sqrt(eb$values)),eb$vectors)
          vec[, k:k1] <- t(tcrossprod(rB,vec[, k:k1]))
      }
      else {
          vec <- t(t(vec)/sqrt(ed$val[1:k]))
      }
      d <- crossprod(vec,R%*%p)
      d <- sum(d^2)
      attr(d, "rank") <- rank
      d
}
  
  lf <- length(coef(object))
  F  <- object$F[1:lf,1:lf]
  Vr <- object$Vb[1:lf,1:lf] 
  table1.1 <- NULL 
  l.sp <- object$l.sp 

  SE <- sqrt( diag(object$Vb) )       
  n  <- object$n 

  estimate1 <- coef(object)[1:(object$gam.fit$nsdf)]
  se1       <- SE[1:(object$gam.fit$nsdf)]
  ratio1    <- estimate1/se1
  pv1       <- 2*pnorm(abs(ratio1), lower.tail = FALSE)
  table1    <- cbind(estimate1,se1,ratio1,pv1)

  dimnames(table1)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")


  if(l.sp!=0){
  	pTerms.df1 <- pTerms.chi.sq1 <- pTerms.pv1 <- edf1 <- NA
		for(k in 1:l.sp){
			ind <- (object$gam.fit$smooth[[k]]$first.para):(object$gam.fit$smooth[[k]]$last.para)
			edf1[k] <- sum(diag(F)[ind])
			names(edf1)[k] <- object$gam.fit$smooth[[k]]$label 
			b  <- coef(object)[ind]
			V  <- Vr[ind,ind]
			Xt <- object$X[, 1:length(ind)+object$gam.fit$nsdf]
			pTerms.df1[k] <- min(ncol(Xt), edf1[k])
			pTerms.chi.sq1[k] <- Tp <- testStat(b, Xt, V, pTerms.df1[k])
			pTerms.df1[k] <- attr(Tp, "rank")
                        pTerms.pv1[k] <- pchisq(pTerms.chi.sq1[k], df = pTerms.df1[k], lower.tail = FALSE)  
                                        }
  	table1.1 <- cbind(edf1, pTerms.df1, pTerms.chi.sq1, pTerms.pv1)
        dimnames(table1.1)[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
                      }


                      res <- list(tableP=table1, 
                                  tableNP=table1.1, 
                                  n=n, tau=object$tau, 
                                  formula=object$gam.fit$formula,
                                  l.sp=l.sp, 
                                  t.edf=object$t.edf)

                      class(res) <- "summary.bgeva"

                                      
res

}



