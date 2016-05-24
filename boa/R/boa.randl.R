"boa.randl" <-
function(link, q, error, prob, delta)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   n <- nrow(link)
   phi <- qnorm(0.5 * (1 + prob))
   n.min <- ceiling(q * (1 - q) * (phi / error)^2)
   if(n.min <= n) {
      pnames <- boa.pnames(link)
      for(i in pnames) {
         dichot <- ifelse(link[, i] <= quantile(link[, i], probs = q), 1, 0)
         thin <- 0
         bic <- 1
         while(bic >= 0) {
            thin <- thin + 1
            test <- dichot[seq(1, n, by = thin)]
            n.test <- length(test)
            temp <- test[1:(n.test - 2)] + 2 * test[2:(n.test - 1)] +
                    4 * test[3:(n.test)]
            tran.test <- array(tabulate(temp + 1, nbins = 8), dim=c(2, 2, 2))
            g2 <- 0
            for(i1 in c(1, 2)) {
               for(i2 in c(1, 2)) {
                  for(i3 in c(1, 2)) {
                     if(tran.test[i1, i2, i3] > 0) {
                        fitted <- log(sum(tran.test[i1, i2, ])) +
                                  log(sum(tran.test[ , i2, i3])) -
                                  log(sum(tran.test[ , i2, ]))
                        g2 <- g2 + 2 * tran.test[i1, i2, i3] *
                                   (log(tran.test[i1, i2, i3]) - fitted)
                     }
                  }
               }
            }
            bic <- g2 - 2 * log(n.test - 2)
         }
         tran.final <- tabulate(test[1:(n.test - 1)] + 2 * test[2:n.test] + 1,
                                nbins = 4)
         alpha <- tran.final[3] / (tran.final[1] + tran.final[3])
         beta <- tran.final[2] / (tran.final[2] + tran.final[4])
         burnin <- ceiling(log(delta * (alpha + beta) / max(alpha, beta)) /
                           log(abs(1 - alpha - beta))) * thin
         keep <- ceiling((2 - alpha - beta) * alpha * beta * phi^2 /
                         (error^2 * (alpha + beta)^3)) * thin
         total <- burnin + keep
         result <- rbind(result, c(thin, burnin, total, n.min, total / n.min))
      }
      dimnames(result) <- list(pnames, c("Thin", "Burn-in", "Total",
                                         "Lower Bound", "Dependence Factor"))
   } else {
      result <- n.min
      names(result) <- "Lower Bound"
   }

   return(result)
}
