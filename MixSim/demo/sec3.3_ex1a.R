### Example 1 of Section 3.3 plot (a)
set.seed(1234)
Q <- MixSim(MaxOmega = 0.20, BarOmega = 0.05, K = 5, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

