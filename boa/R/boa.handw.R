"boa.handw" <-
function(link, error, alpha)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   pnames <- boa.pnames(link)
   stest <- keep <- I <- htest <- xbar <- halfwidth <-
            structure(rep(NA, length(pnames)), names = pnames)
   iter <- unique(boa.iter(link))
   n <- length(iter)
   n.min <- round(0.5 * n)
   idx.drop <- -1 * 1:max(round(0.10 * n), 1)
   S0 <- boa.gewekePwr(boa.getiter(link, iter[(n - n.min + 1):n]))
   q.upper <- qnorm(1 - alpha / 2)
   for(i in pnames) {
      parm <- boa.getparms(link, i)
      piter <- iter
      keep[i] <- n
      stest[i] <- "failed"
      htest[i] <- "failed"
      while((keep[i] >= n.min) && (stest[i] == "failed")) {
         n.parm <- nrow(parm)
         xbar[i] <- mean(parm)
         halfwidth[i] <- q.upper * sqrt(boa.gewekePwr(parm) / n.parm)
         if(abs(halfwidth[i] / xbar[i]) <= error)  htest[i] <- "passed"
         B <- cumsum(parm) - xbar[i] * 1:n.parm
         Bsq <- (B * B) / (n.parm * S0[i])
         I[i] <- sum(Bsq) / n.parm
         if(I[i] < 0.46) {
            stest[i] <- "passed"
         } else {
            piter <- piter[idx.drop]
            keep[i] <- length(piter)
            parm <- boa.getiter(parm, piter)
         }
      }
   }
   result <- data.frame(stest, keep, n - keep, I, htest, xbar, halfwidth)
   names(result) <- c("Stationarity Test", "Keep", "Discard", "C-von-M",
                      "Halfwidth Test", "Mean", "Halfwidth")

   return(result)
}
