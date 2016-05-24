`plot.adonis.pertables` <-
function(x, ...) {
n <- dim(x$simulation$p.quant)[2]
layout(matrix(c(1:(n*2)), n, 2))
for(i in 1:n){
par(mar=c(4, 4, 3, 1))
z <- x$simulation$R2[i, ]
diff.z <- diff(c(min(z), max(z)))
hist(z, breaks = seq(min(z), max(z), by = diff.z/30), main = "", xlab = "R2 (%)", col = "lightblue", border = "darkred", xlim = c(min(z)-diff.z*0.25, max(z)+diff.z*0.25))
mtext(rownames(x$simulation$R2)[i], adj = 0, cex = 1.2, line = 1)
par(new = TRUE)
plot(density(z), col = "red", lty = 3, lwd = 2, main = "", xlab = "", ylab = "", axes = FALSE, xlim = c(min(z)-diff.z*0.25, max(z)+diff.z*0.25))
abline(v = median(z), lwd = 2, col = "red")
}
for(i in 1:n){
par(mar=c(4, 4, 3, 1))
z <- x$simulation$pvalue[i, ]
diff.z <- diff(c(min(z), max(z)))
a<- hist(z, plot = FALSE, breaks = seq(min(z), max(z), by = diff.z/30))$mids
hist(z, breaks = seq(min(z), max(z), by = diff.z/30), main = "", xlab = "p-value", col = ifelse(a < 0.05, "lightblue", "white"), border = "darkred", xlim = c(min(z)-diff.z*0.25, max(z)+diff.z*0.25))
mtext(rownames(x$R2)[i], adj = 0, cex = 1.2, line = 1)
par(new = TRUE)
plot(density(z), col = "red", lty = 3, lwd = 2, main = "", xlab = "", ylab = "", axes = FALSE, xlim = c(min(z)-diff.z*0.25, max(z)+diff.z*0.25))#, add=TRUE)
abline(v = median(z), lwd = 2, col = "red")
}
}

