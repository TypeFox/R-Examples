contord<-function (marginal, Sigma, support = list(), Spearman = FALSE) 
{
    if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
    if(!isSymmetric(Sigma) | min(eigen(Sigma)$values)<0 | !all(diag(Sigma)==1)) stop("Correlation matrix not valid!")
    len <- length(support)
    k <- length(marginal)
    kj <- numeric(k)
    Sigmaord <- diag(k)
    for (i in 1:k) {
        kj[i] <- length(marginal[[i]]) + 1
        if (len == 0) {
                      support[[i]] <- 1:kj[i]
                      }
        if (Spearman) {
	              s1 <- c(marginal[[i]], 1)
	              s2 <- c(0,marginal[[i]])
                      support[[i]] <- (s1+s2)/2
                      }
                    }
    L <- vector("list", k)
    for (i in 1:k) {
                   L[[i]] <- qnorm(marginal[[i]])
                   L[[i]] <- c(-Inf, L[[i]], +Inf)
                   }
    for (q in 1:(k - 1)) {
                         for (r in (q + 1):k) {
                         pij <- matrix(0, kj[q], kj[r])
                         for (i in 1:kj[q]) {
                         for (j in 1:kj[r]) {
                                            low <- rep(-Inf, k)
					    upp <- rep(Inf, k)
                                            low[q] <- L[[q]][i]
                                            low[r] <- L[[r]][j]
                                            upp[q] <- L[[q]][i + 1]
                                            upp[r] <- L[[r]][j + 1]
                                            pij[i, j] <- pmvnorm(low, upp, rep(0, k), corr = Sigma)
                                            low <- rep(-Inf, k)
                                            upp <- rep(Inf, k)
                                            }
                                            }
                         my <- sum(apply(pij, 2, sum) * support[[r]])
                         sigmay <- sqrt(sum(apply(pij, 2, sum) * support[[r]]^2) - my^2)
                         mx <- sum(apply(pij, 1, sum) * support[[q]])
                         sigmax <- sqrt(sum(apply(pij, 1, sum) * support[[q]]^2) - mx^2)
                         mij <- support[[q]] %*% t(support[[r]])
                         muij <- sum(mij * pij)
                         covxy <- muij - mx * my
                         corxy <- covxy/(sigmax * sigmay)
                         Sigmaord[q, r] <- corxy
                          }
    }
    as.matrix(forceSymmetric(Sigmaord))
}
