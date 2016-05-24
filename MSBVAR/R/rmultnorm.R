"rmultnorm" <-
function(n, mu, vmat, tol = 1e-10)
       {
               p <- ncol(vmat)
               if(length(mu)!=p)
                       stop(paste("mu vector is the wrong length:",length(mu)))
#               if(max(abs(vmat - t(vmat))) > tol)
#                       stop("vmat not symmetric")
               vs <- La.svd(vmat)
               vsqrt <- t(t(vs$vt) %*% (t(vs$u) * sqrt(vs$d)))
               ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
               ans <- sweep(ans, 2, mu, "+")
               dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])
               ans
       }

