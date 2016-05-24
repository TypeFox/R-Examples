### Example 2 of Section 3.3 plot (c)
set.seed(1238)
Q <- MixSim(MaxOmega = 0.1, K = 3, p = 2, int = c(0.2, 1))
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, lambda = c(0.1, 10))
colors <- c("red", "green", "blue")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

