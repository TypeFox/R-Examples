# A = pJ x pJ Fisher information matrix (based on the 2nd derivatives of the log-likelihood with respect to B)
# W = pJ x [pJ]^2 matrix of 3rd derivatives
# V = pJ x [pJ]^3 matrix of 4th derivatives

getV <- function(x, wt, B) {

   p <- nrow(B) # p covariates
   J <- ncol(B) # J categories

   f1 <- texp(x %*% B)
   f2 <- f1/(1 + apply(f1, 1, sum))

   A <- matrix(data = 0, nrow = p*J, ncol = p*J)
   W <- matrix(data = 0, nrow = p*J, ncol = (p*J)^2)
   V <- matrix(data = 0, nrow = p*J, ncol = (p*J)^3)

   ij <-1
   while (ij <= J) {
      ik <- 1
      while (ik <= J) {
         if (ij == ik) {
            s <- f2[,ij] * (1 - f2[,ik]); cc1 <- 1
         } else {
            s <- f2[,ij] * f2[,ik]; cc1 <- -1
         }
         sx <- sqrt(s) * x * sqrt(wt)
         A <- A + cc1 * (crossprod(sx) %x% (cbind(diag(J)[,ij]) %*% rbind(diag(J)[ik,])))
         il <- 1
         while (il <= J) {
            if (ij != ik & ij != il & ik != il) {
               q <- 2 * f2[,ij] * f2[,ik] * f2[,il]; cc2 <- 1
            }
            if (ij == ik & ij == il) {
               q <- f2[,ij] * (1 - f2[,ij]) * (1 - 2 * f2[,ij]); cc2 <- 1
            } else {
               if (ij == ik) {
                  q <- f2[,ij] * f2[,il] * (1 - 2 * f2[,ij]); cc2 <- -1
               }
               if (ij == il) {
                  q <- f2[,ij] * f2[,ik] * (1 - 2 * f2[,ij]); cc2 <- -1
               }
               if (ik == il) {
                  q <- f2[,ij] * f2[,ik] * (1 - 2 * f2[,ik]); cc2 <- -1
               }
            }
            qxab <- t(wt * q * x) %*% (x %*~% x)
            W <- W + (cc2 * qxab %*% (diag(p) %x% rbind(diag(J)[ik,]) %x% diag(p)) %x% (cbind(diag(J)[,ij]) %*% rbind(diag(J)[il,])))
            im <- 1
            while (im <= J) {
               if (ij == ik & ij == il & ij == im) {
                  pp <- f2[,ij] * (1 - f2[,ij]); ss <- pp * (1 - 6 * pp); cc3 <- 1
               } else {
                  if (ij == ik & ij == il) {
                     ss <- f2[,ij] * f2[,im] * (1 - 6 * f2[,ij] * (1 - f2[,ij])); cc3 <- -1
                  }
                  if (ij == ik & ij == im) {
                     ss <- f2[,ij] * f2[,il] * (1 - 6 * f2[,ij] * (1 - f2[,ij])); cc3 <- -1
                  }
                  if (ij == il & ij == im) {
                     ss <- f2[,ij] * f2[,ik] * (1 - 6 * f2[,ij] * (1 - f2[,ij])); cc3 <- -1
                  }
                  if (ik == il & ik == im) {
                     ss <- f2[,ij] * f2[,ik] * (1 - 6 * f2[,ik] * (1 - f2[,ik])); cc3 <- -1
                  }
                  if (ij == ik & il == im) {
                     pp <- f2[,ij] * f2[,il]; ss <- pp * (1 - 2 * (f2[,ij] + f2[,il]) + 6 * pp); cc3 <- -1
                  }
                  if (ij == il & ik == im) {
                     pp <- f2[,ij] * f2[,ik]; ss <- pp * (1 - 2 * (f2[,ij] + f2[,ik]) + 6 * pp); cc3 <- -1
                  }
                  if (ij == im & ik == il) {
                     pp <- f2[,ij] * f2[,ik]; ss <- pp * (1 - 2 * (f2[,ij] + f2[,ik]) + 6 * pp); cc3 <- -1
                  }
                  if (ij == ik & ij != il & ij != im & ik != il & ik != im & il != im) {
                     ss <- f2[,ij] * f2[,il] * f2[,im] * (1 - 3 * f2[,ij]); cc3 <- 2
                  }
                  if (ij != ik & ij == il & ij != im & ik != il & ik != im & il != im) {
                     ss <- f2[,ij] * f2[,ik] * f2[,im] * (1 - 3 * f2[,ij]); cc3 <- 2
                  }
                  if (ij != ik & ij != il & ij != im & ik == il & ik != im & il != im) {
                     ss <- f2[,ij] * f2[,ik] * f2[,im] * (1 - 3 * f2[,ik]); cc3 <- 2
                  }
                  if (ij != ik & ij != il & ij != im & ik != il & ik == im & il != im) {
                     ss <- f2[,ij] * f2[,ik] * f2[,il] * (1 - 3 * f2[,ik]); cc3 <- 2
                  }
                  if (ij != ik & ij != il & ij == im & ik != il & ik != im & il != im) {
                     ss <- f2[,ij] * f2[,ik] * f2[,il] * (1 - 3 * f2[,ij]); cc3 <- 2
                  }
                  if (ij != ik & ij != il & ij != im & ik != il & ik != im & il == im) {
                     ss <- f2[,ij] * f2[,ik] * f2[,il] * (1 - 3 * f2[,il]); cc3 <- 2
                  }
               }
               if (ij != ik & ij !=il & ij != im & ik != il & ik != im & il != im) {
                  ss <- f2[,ij] * f2[,ik] * f2[,il] * f2[,im]; cc3 <- -6
               }
               ssxab <- t(wt * ss * x) %*% (x %*~% x %*~% x)
               V <- V + cc3 * (ssxab %*% (diag(p) %x% rbind(diag(J)[ik,]) %x% diag(p) %x% rbind(diag(J)[il,]) %x% diag(p)) %x% (cbind(diag(J)[,ij]) %*% rbind(diag(J)[im,])))
               im <- im + 1
            }
            il <- il + 1
         }
         ik <- ik + 1
      }
      ij <- ij + 1
   }

   list(A = A, W = W, V = V)

}
