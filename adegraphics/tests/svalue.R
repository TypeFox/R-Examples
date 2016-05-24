library(adegraphics)
pdf("svalue.pdf")

## ex1
xy <- cbind.data.frame(x = runif(50, -1, 1), y = runif(50, 0, 2))
z <- rnorm(50)
z <- sapply(z, function(X) max(X, -3))
z <- sapply(z, function(X) min(X, 3))
val1 <- s.value(xy, z, method = "size", symbol = "square", plot = F)
val2 <- s.value(xy, z, method = "color", symbol = "square", plot = F)
val3 <- s.value(xy, z, method = "size", symbol = "circle", plot = F)
val4 <- s.value(xy, z, method = "color", symbol = "circle", plot = F)
g1 <- ADEgS(c(val1, val2, val3, val4), positions = layout2position(matrix(c(1, 2, 3, 4), 2, 2)), add = matrix(0, ncol = 4, nrow = 4))
g2 <- s.value(xy, z, method = "color", symbol = "square", breaks = c(-3, -1, -0.5, 0, 0.5, 1, 3))
g3 <- s.value(xy, z, method = "color", col = colorRampPalette(c("yellow", "blue"))(6))
g4 <- s.value(xy, z, method = "size", symbol = "circle", paxes.draw = FALSE)

## ex2
xx <- runif(100) * 100
yy <- 1:100 
zz <- 1:100
breaks <- c(0, 25, 50, 75, 100)
g5 <- s.value(data.frame(xx, yy), zz, breaks = breaks, method = "color", paxes.draw = TRUE,  porigin.include = FALSE)
g6 <- s.value(data.frame(xx, yy), cbind(zz, rev(zz)), breaks = breaks, method = "color", col = c("blue", "red", "green", "yellow"), paxes.draw = TRUE)
g7 <- s.value(data.frame(xx, yy), cbind(zz, rev(zz)), nclass = c(2, 6), method = "color", col = c("blue", "red", "pink", "green", "yellow"), paxes.draw = TRUE)

## ex3
data(rpjdl, package = "ade4")
fau.coa <- ade4::dudi.coa(rpjdl$fau, scan = FALSE, nf = 3)
val5 <- s.value(fau.coa$li, fau.coa$li[,3], plot = FALSE)
val6 <- s.value(fau.coa$li, fau.coa$li[, 3], center = 0, method = "size", symbol = "circle", col = c("yellow", "red"), plot = FALSE)
g8 <- ADEgS(c(val5, val6), positions = layout2position(matrix(c(1, 2), 1, 2)), add = matrix(0, ncol = 2, nrow = 2))

## ex3
data(doubs, package = "ade4")
g9 <- s.value(doubs$xy, doubs$env[, 1:2])
g10 <- s.value(doubs$xy, doubs$env)
