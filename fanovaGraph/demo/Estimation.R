#####################################################################
### Comparisons of estimation methods for total interaction indices #
#####################################################################

require(fanovaGraph)

### function definition
domain <- c(0, 1)
d <- 6
a <- c(0, 0, 0, 2/5, 2/5, 5)
fun <- function(X) {
    y <- 1
    for (j in 1:d) {
        y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
    }
    y
}

true <- c(0.20470157, 0.20470157, 0.119012525, 0.119012525, 
    0.007511985, 0.20470157, 0.119012525, 0.119012525, 0.007511985, 0.119012525, 
    0.119012525, 0.007511985, 0.069193319, 0.004367432, 0.004367432)

### estimations

N.run <- 10
N.eval <- 20000
Int1 <- matrix(, N.run, choose(d, 2))
Int2 <- matrix(, N.run, choose(d, 2))
Int3 <- matrix(, N.run, choose(d, 2))
Int4 <- matrix(, N.run, choose(d, 2))

for (i in 1:N.run) {
    print(paste("i=", i))
    Int1[i, ] <- estimateGraph(fun, d = d, n.tot = N.eval, method = "LiuOwen", 
        q.arg = list(min = domain[1], max = domain[2]))$tii[,1]
    Int2[i, ] <- estimateGraph(fun, d = d, n.tot = N.eval, method = "FixFast", 
        q.arg = list(min = domain[1], max = domain[2]))$tii[,1]
    Int3[i, ] <- estimateGraph(fun, d = d, n.tot = N.eval, method = "RBD", 
        q.arg = list(min = domain[1], max = domain[2]))$tii[,1]
    Int4[i, ] <- estimateGraph(fun, d = d, n.tot = N.eval, method = "PickFreeze", 
        q.arg = list(min = domain[1], max = domain[2]))$tii[,1]
}

### boxplots

plot(0, type = "n", ylim = c(min(Int1, Int2, Int3, Int4) * 1.1, 
    max(Int1, Int2, Int3, Int4) * 1.1), xlim = c(1, choose(d, 2)), xaxt = "n", 
    xlab = "interaction", ylab = "total interaction index estimation")
boxplot(Int1, add = TRUE, at = 1:choose(d, 2) - 0.3, 
    boxwex = 0.1, xaxt = "n", pch = 3, cex = 0.5)
boxplot(Int2, add = TRUE, at = 1:choose(d, 2) - 0.1, col = 2, 
    boxwex = 0.1, xaxt = "n", pch = 3, cex = 0.5)
boxplot(Int3, add = TRUE, at = 1:choose(d, 2) + 0.1, col = 3, 
    boxwex = 0.1, xaxt = "n", pch = 3, cex = 0.5)
boxplot(Int4, add = TRUE, at = 1:choose(d, 2) + 0.3, col = 4, 
    boxwex = 0.1, xaxt = "n", pch = 3, cex = 0.5)
points(1:choose(d, 2), true, cex = 1, pch = 4, col=1)
abline(h = 0, v = 1:(choose(d, 2) - 1) + 0.5, lty = 3)
axis(1, at = 1:choose(d, 2), labels = paste(combn(d, 2)[1, 
    ], combn(d, 2)[2, ], sep = ""))
legend("topright", legend = c("LiuOwen", "FixFast", "RBD",
    "PickFreeze", "true value"), pch = c(22, 22, 22, 22, 4), col = 1, cex = 0.6, 
    pt.bg = c(0, 2,3,4)) 