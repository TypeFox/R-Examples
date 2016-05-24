CGaPred <-
function(M, xf = "median", confidence = 0.95) {
  K <- M$K
  p <- M$p
  tao <- M$tao
  covar <- M$covar
  MS <- M$summary
  prob <- (1 - confidence) / 2
  SUM <- CLambdaSumm(M, confidence)
  SUM.h <- SUM$SUM.h
  h.0 <- SUM$SUM.h[, 2]
  S.0 <- SUM$SUM.S[, 2]
  H.0 <- SUM$SUM.H[, 2]
  b <- 0
  theta.summary <- matrix(0, ncol = 5, nrow = p)
  b <- 0
  if(length(MS[, 1]) == 3 * K - 1 + p){
    b <- 1
  }
  for(i in 1:p) {
    theta.summary[i, 1] <- mean(MS[3 * K - 2 + b + i, ])
    theta.summary[i, 2] <- quantile(MS[3 * K - 2 + b + i, ], probs = prob)
    theta.summary[i, 3] <- quantile(MS[3 * K - 2 + b + i, ], probs = 1 - prob)
    theta.summary[i, 4] <- median(MS[3 * K - 2 + b + i, ])
    theta.summary[i, 5] <- sd(MS[3 * K - 2 + b + i, ])
  }
  colnames(theta.summary) <- c("mean", prob, 1- prob, "median", "sd")
  theta <- theta.summary[, 1]
  if (class(xf) == "character") {
    xf <- rep(0, p)
    for (i in 1:p) {
      xf[i] <- median(covar[, i])  
    }
  }
  h.xf <- h.0 * exp(theta %*% xf)
  S.xf <- exp(- H.0 * exp(theta %*% xf))
  plot(c(0, max(tao)), c(0, max(h.xf)), "n", xlab = "times", ylab = "", 
       main = "Hazard distribution estimate")
  for(i in 1:K) {
    segments(x0 = tao[i], y0 = SUM.h[i, 2], x1 = tao[i + 1], 
             y1 = SUM.h[i, 2], lty = 1, lwd = 2.5)
    segments(x0 = tao[i], y0 = h.xf[i], x1 = tao[i + 1], 
             y1 = h.xf[i], lty = 1, lwd = 2.5, col="red")
  }
  legend("bottomright", c("Baseline hazard", "Estimate for x_f"), 
         lty = c(1, 1), lwd = c(2.5, 2.5), col = c(1, "red"), bty = "n",
         cex = 0.8)
  out <- list(theta.summary = theta.summary, h.xf = h.xf, S.xf = S.xf)
  return(out)
}
