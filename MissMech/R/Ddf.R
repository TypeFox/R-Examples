Ddf <- function(data, mu, sig) {
 y <- data
 n <- nrow(y)
 p <- ncol(y)
 ns <- p * (p + 1) / 2
 nparam <- ns + p
 ddss <- matrix(0, ns, ns)
 ddmm <- matrix(0, p, p)
 ddsm <- matrix(0, ns, p)
 for (i in 1:n) {
     obs <- which(!is.na(y[i, ]))
     lo <- length(obs)
     tmp <- cbind(rep(obs, 1, each = lo), rep(obs, lo))
     tmp <- matrix(tmp[tmp[, 1] >= tmp[, 2], ], ncol = 2)
     lolo = lo*(lo+1)/2
     subsig <- sig[obs, obs]
     submu <- mu[obs, ]
     temp <- matrix(y[i, obs] - submu, nrow = 1)
     a <- solve(subsig)
     b <- a %*% (2 * t(temp) %*% temp - subsig) * a
     d <- temp %*% a
     ddimm <- 2 * a
#==================== DD(mu, mu)
     ddmm[obs, obs] <- ddmm[obs, obs] + ddimm 
#==================== DD(sig, mu)
     rcnt <- 0
     ddism <- matrix(0, lolo, lo)
     for (k in 1:lo) {
         for (l in 1:k) {
             rcnt <- rcnt + 1
             ccnt <- 0
             for (kk in 1:lo) {
                 ccnt <- ccnt + 1
                 ddism[rcnt, ccnt] <- 2 * (1 - 0.5 * (k == l)) * 
                                      (a[kk, l] %*% d[k] + a[kk, k] %*% d[l])
             }
         }
     }
     for (k in 1:lolo) {
         par1 <- tmp[k, 1] * (tmp[k, 1] -1) / 2 + tmp[k, 2]
         for (j in 1:lo) {
             ddsm[par1, obs[j]] <- ddsm[par1, obs[j]] + ddism[k, j]  
         } 
     }
#------------------------ Test part
     ssi <- matrix(0, lolo, lolo)
     for (m in 1:lolo) {   
         u <- which(obs == tmp[m, 1])
         v <- which(obs == tmp[m, 2])
         for (q in 1:m) {
         k <- which(obs == tmp[q, 1])
         l <- which(obs == tmp[q, 2])
         ssi[m, q] <- (b[v, k] * a[l, u] + b[v, l] * a[k, u] + b[u, k] * a [l, v] +
                      b[u, l] * a[k, v]) * (1 - 0.5 * (u == v)) * (1 - 0.5 * (k == l)) 
         }
     }   
     for (k in 1:lolo) {
         par1 <- tmp[k, 1] * (tmp[k, 1] -1) / 2 + tmp[k, 2]
         for (l in 1:k) {
         par2 <- tmp[l, 1] * (tmp[l, 1] - 1) / 2 + tmp[l, 2]
         ddss[par1, par2] <- ddss[par1, par2] + ssi[k, l]
         ddss[par2, par1] <- ddss[par1, par2]
         }
     }
 }
 dd <- -1 * rbind(cbind(ddmm, t(ddsm)), cbind(ddsm, ddss)) / 2
 se <- -solve(dd)
 list(dd = dd, se = se)
}
