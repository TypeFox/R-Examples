### Example 4 of Section 3.3 plot (d)
set.seed(1236)
Q <- MixSim(MaxOmega = 0.1, K = 4, p = 1)
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 1)
colors <- c("red", "green", "blue", "brown")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

