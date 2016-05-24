PlotTheta <- function(M, i = 1, plot.all = TRUE, quantiles = TRUE) {
  K <- M$K
  p <- M$p
  MS <- M$summary
  b <- 0
  c <- 1
  d <- p
  quant <- matrix(0, ncol = 5, nrow = 2)
  m_sd <- matrix(0, ncol = 2, nrow = 2)
  if(length(MS[, 1]) == 3 * K - 1 + p){
    b <- 1
  }
  if (plot.all == FALSE) {
    c <- d <- i
  }
  for(s in c:d) {
    X <- MS[3 * K - 2 + b + s, ]
    hist(X, prob = TRUE, main = paste("Histogram and density for theta_", s, 
                                      sep = ""), 
         xlab = paste("theta_", s, sep = ""), col = "skyblue")
    lines(density(X), lwd = 2)
    legend(x = "topright", legend = c("Histogram", "Density"), lty = c(0, 0), 
           col = c("skyblue", "white"), bty = "n", cex = 0.8, 
           fill = c("skyblue", 1))
    if(s < d) {
      par(mfrow = c(1, 1), ask = TRUE)
    }
    if (quantiles == TRUE) {
      quant[s, ] <- quantile(X, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
      m_sd[s, 1] <- mean(X)
      m_sd[s, 2] <- sd(X)
    }
  }
  quant <- cbind(m_sd, quant)
  quant <- as.data.frame(quant)
  names(quant) <- c("mean", "sd", names(quantile(1:100, 
                                     probs = c(0.025,0.05,0.5,0.95,0.975))))
  for (s in c:d) {
    row.names(quant)[s] <- paste("theta_", s, sep = "")
  } 
  par(mfrow = c(1, 1), ask = FALSE)
  return(quant)
}

