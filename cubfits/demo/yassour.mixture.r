suppressMessages(library(cubfits, quietly = TRUE))
library(EMCluster)
set.seed(1234)

### Get individual of phi.Obs.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs <- GM / sum(GM) * 15000

### Run EM optimization.
X <- log(as.matrix(phi.Obs))
K.max <- 6
ret <- list()
for(i.k in 1:K.max){
  ret[[i.k]] <- rand.EM(X, nclass = i.k)
  tmp <- summary(ret[[i.k]])
  cat("K = ", i.k, ", logL = ", tmp$logL,
      ", AIC = ", tmp$AIC, ", BIC = ", tmp$BIC, "\n", sep = "")
}

### GM averaged
x <- seq(min(X), max(X), length = 100)
d.em <- list()
for(i.k in 1:K.max){
  d.em[[i.k]] <- sapply(x, dmixmvn, ret[[i.k]])
}
ylim <- c(0, 0.5)

hist(X, nclass = 50, freq = FALSE, ylim = ylim,
     xlab = "Expression (log)", ylab = "Density",
     main = paste("Yassour (N = ", length(X), ", averaged)",
                  sep = ""))
for(i.k in 1:K.max){
  lines(x = x, y = d.em[[i.k]], col = i.k, lty = i.k, lwd = 1.5)
}
legend(min(x) + 0.5, ylim[2] * 0.9,
       paste("K = ", 1:K.max, sep = ""),
       col = 1:K.max, lty = 1:K.max, lwd = 1.5)

