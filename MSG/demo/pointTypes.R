par(mar = c(0.5, 0.5, 0.5, 0.5))
ipch = 0:33; np = length(ipch); k = floor(sqrt(np))
dd = c(-1, 1)/2; rx = dd + range(ix <- ipch%/%k)
ry = dd + range(iy <- 3 + (k - 1) - ipch%%k)
pch = as.list(ipch)
pch[26 + 1:8] = as.list(c("*", ".", "o", "O", "0",
    "+", "-", "|"))
plot(rx, ry, type = "n", axes = FALSE, ann=FALSE)
abline(v = ix, h = iy, col = "gray", lty = "dotted")
for (i in 1:np) {
   points(ix[i], iy[i], pch = pch[[i]], bg = "yellow", cex = 3)
   text(ix[i] - 0.3, iy[i], pch[[i]], col = "brown", cex = 1.2)
}
box()
