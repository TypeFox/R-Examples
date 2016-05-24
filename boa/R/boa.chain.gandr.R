"boa.chain.gandr" <-
function(chain, chain.support, alpha, pnames, window, to)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   lnames <- names(chain)
   if(missing(pnames)) {
      pnames <- boa.pnames(chain[[1]])
   } else {
      pnames <- intersect(boa.pnames(chain[[1]]), pnames)
   }
   iter <- boa.iter(chain[[1]])
   if(missing(window))  window <- 0.5
   if(!missing(to))  iter <- iter[iter <= to]
   for(i in lnames[-1]) {
      pnames <- intersect(boa.pnames(chain[[i]]), pnames)
      iter <- intersect(boa.iter(chain[[i]]), iter)
   }
   n <- length(iter)
   keep <- iter[(n - round(window * n) + 1):n]

   m <- length(lnames)
   n <- length(keep)
   p <- length(pnames)
   q.upper <- 1 - alpha / 2

   if(m < 2) {
      cat("Warning: must supply at least two sequences to analyze\n")
   } else if(p == 0) {
      cat("Warning: no common parameters to analyze\n")
   } else if(n < (p + 1)) {
      cat("Warning: must supply at least", 2 * p + 1,
          "common iterations to analyze\n")
   } else {
      psi <- matrix(NA, nrow = n, ncol = p, dimnames = list(NULL, pnames))
      xbar.within <- sxx.within <- matrix(NA, nrow = m, ncol = p,
                                          dimnames = list(lnames, pnames))
      W <- 0
      for(i in lnames) {
         rows <- is.element(boa.iter(chain[[i]]), keep)
         for(j in pnames) {
            psi[, j] <- boa.transform(chain[[i]][rows, j],
                                      chain.support[[i]][, j])
         }
         xbar.within[i, ] <- colMeans(psi)
         sxx <- var(psi)
         sxx.within[i, ] <- diag(sxx)
         W <- W + sxx
      }
      W <- W / m
      B <- n * var(xbar.within)

      w <- diag(W)
      b <- diag(B)
      xbar <- colMeans(xbar.within)
      var.w <- colVars(sxx.within) / m
      df.w <- 2 * w^2 / var.w
      var.b <- 2 * b^2 / (m - 1)
      var.wb <- NULL
      for(j in pnames) {
         var.wb <- c(var.wb,
                     var(sxx.within[, j], xbar.within[, j]^2) -
                     2 * xbar[j] * var(sxx.within[, j], xbar.within[, j]))
      }
      var.wb <- (n / m) * var.wb
      v <- ((n - 1) * w + (1 + 1 / m) * b) / n
      var.v <- ((n - 1)^2 * var.w + (1 + 1 / m)^2 * var.b + 2 * (n - 1) *
                (1 + 1 / m) * var.wb) / n^2
      df.v <- 2 * v^2 / var.v

      psrf <- sqrt(v / w)
      names(psrf) <- pnames
      csrf <- sqrt((df.v + 3) / (df.v + 1) * cbind(v / w, (1 - 1 / n) +
                   qf(q.upper, m - 1, df.w) * (1 + 1 / m) * b / (n * w)))
      dimnames(csrf) <- list(pnames, c("Estimate", q.upper))
      mpsrf <- sqrt((1 - 1 / n) + (1 + 1 / p) * eigen(qr.solve(W, B / n),
                    symmetric = FALSE, only.values = TRUE)$values[1])
      result <- list(psrf = psrf, csrf = csrf, mpsrf = mpsrf,
                     window = range(keep))
   }
   if(is.null(result)) {
      result <- list(psrf = structure(rep(NA, p), names = pnames),
                     csrf = matrix(NA, nrow = p, ncol = 2,
                            dimnames = list(pnames, c("Estimate", q.upper))),
                     mpsrf = NA, window = c(NA, NA))
   }

   return(result)
}
