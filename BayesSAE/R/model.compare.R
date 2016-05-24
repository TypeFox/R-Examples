model.compare <-
function(result, subset, poest = "mean"){
     innov <- result$innov
     if (innov != "normal" && innov != "t")
          stop("Sampling model must be with normal or t innovations") 
     if (is.null(result$theta))
          stop("MCMC draws of theta's are needed")
     if (poest != "mean" && poest != "median")
          stop("only posterior means or medians can be used as the point estimator for theta's at present")
     theta <- result$theta[,subset]
     m <- nrow(theta)
     n <- ncol(theta)
     Y <- result$Y
     if (length(Y) != m)
          stop("length of response vector and number of rows of theta should be the same")
     if (innov == "normal"){
          Z <- result$Z
          if (length(Z) != m)
               stop("length of vardir and number of rows of theta should be the same")
          D_avg <- -2 * mean(apply(dnorm(Y, theta, sqrt(Z), log = TRUE), 2, sum))
          if (poest == "mean")
               D_theta.hat <- dnorm(Y, rowMeans(theta), sqrt(Z), log = TRUE)
          else
               D_theta.hat <- dnorm(Y, apply(theta, 1, mean), sqrt(Z), log = TRUE) 
          D_theta.hat <- -2 * sum(D_theta.hat)  
     }
     else{
          sig2 <- result$sig2[,subset]
          if (nrow(sig2) != m)
               stop("number of rows of sigma2 should be the same as that of theta")
          if (ncol(sig2) != n)
               stop("number of draws of sigma2 should be the same as that of theta")
          D_avg <- -2 * mean(apply(dnorm(Y, theta, sqrt(sig2), log = TRUE), 2, sum))
          if (poest == "mean")
               D_theta.hat <- dnorm(Y, rowMeans(theta), sqrt(rowMeans(sig2)), log = TRUE)
          else
               D_theta.hat <- dnorm(Y, apply(theta, 1, mean), sqrt(apply(sig2, 1, mean)), log = TRUE) 
          D_theta.hat <- -2 * sum(D_theta.hat) 
     }
     criter <- list(D_avg = D_avg, D_theta.hat = D_theta.hat, DIC = 2 * D_avg - D_theta.hat)
}
