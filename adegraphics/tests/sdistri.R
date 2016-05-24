library(adegraphics)
pdf("sdistri.pdf")

xy5 <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -1, 1))
w1 <- as.numeric((xy5$x > 0) & (xy5$y > 0))
w2 <- ((xy5$x > 0) & (xy5$y < 0)) * (1 - xy5$y) * xy5$x
w3 <- ((xy5$x < 0) & (xy5$y > 0)) * (1 - xy5$x) * xy5$y
w4 <- ((xy5$x < 0) & (xy5$y < 0)) * xy5$y * xy5$x
distri <- data.frame(a = w1 / sum(w1), b = w2 / sum(w2), c = w3 / sum(w3), d = w4 / sum(w4))
g5 <- s.distri(xy5, distri, plabels.boxes = list(col = "white", alpha = 1), plabels.cex = 2, plabels.col = 1:5)

data(rpjdl, package = "ade4")
xy6 <- ade4::dudi.coa(rpjdl$fau, scan = FALSE)$li + 3
g6 <- s.distri(xy6, rpjdl$fau[, 5], ellipseSize = 1.5, psub = list(text = rpjdl$frlab[5], cex = 2, pos = c(0.2, 0.1)))
g7 <- s.distri(xy6, rpjdl$fau[, 5], ellipseSize = 1.5, psub = list(text = rpjdl$frlab[5], cex = 2, pos = c(0.2, 0.1)), porigin = list(include = FALSE), pellipses.axes.col = "blue")

## test add
g8 <- s.distri(xy6, rpjdl$fau[, 5], ellipseSize = 1.5, psub = list(text = rpjdl$frlab[5], cex = 2, pos = c(0.2, 0.1)), porigin.include = FALSE, pellipses = list(col = "blue"))
g9 <- s.distri(xy6, rpjdl$fau[, 12], ellipseSize = 1.5, psub = list(text = rpjdl$frlab[5], cex = 2, pos = c(0.2, 0.1)), porigin.include = FALSE, pellipses = list(col = "red"), add = TRUE)
show(g9) ## g8 is a superposition, an ADEgS object 

## add
index <- c(1, 5, 8, 20, 21, 23, 26, 33, 36, 44, 47, 49)
col <- colorRampPalette(c("blue", "red", "orange"))(49)
s.distri(xy6, rpjdl$fau[, 1], ellipseSize = 1, starSize = 0, porigin.include = FALSE, pellipses = list(col = col[1], alpha = 0.3))
for(i in index[-1])
   s.distri(xy6, rpjdl$fau[, i], ellipseSize = 1, starSize = 0, porigin.include = FALSE, pellipses = list(col = col[i], alpha = 0.3), add = TRUE)

current <- get("currentadeg", env = adegraphics:::.ADEgEnv)
print(current[[6]])
length(current) == length(index)