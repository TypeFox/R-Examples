library(adegraphics)
pdf("sclass.pdf")

xy0 <- cbind.data.frame(x = runif(20, -1, 1), y = runif(20, -1, 6))
basic <- s.class(xy0, fac = factor(rep(c("A", "B"), le = 20)), chull = 0, star = 0)

xy1 <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -1, 6))
fac1 <- factor(xy1$x > 0) : factor(xy1$y > 0)
g1 <- s.class(xy1, fac = fac1, storeData = F, col = 1:4, pbackground.box = T, pbackground.col = grey(0.85), paxes.draw = T, ell = 0)

## multiaxis
xy2 <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -2, 2), y2 = runif(200, -0.5, .5))
fac2 <- factor(xy2$x > 0) 
g2 <- s.class(xy2, fac = fac2, xax = 1, yax = 2:3, storeData = F, plot = F)
print(g2)

## insertion
print(ADEgS(list(g1, g2), posi = rbind(c(0, 0, 1, 1), c(0.7, 0.5, 1, 0.9))))

## color test
g3 <- s.class(xy1, fac = fac2, psub.text = "Graphic 3", "ppoints.col" = 1:5, plabels.boxes = list(col = "white", alpha = 0.8), pellipses.col = 1:5, col = 1:5)

## test convex hull and parameters
xy4 <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -1, 1))
fac4 <- factor(xy4$x > 0) : factor(xy4$y > 0)
col <- c("black", "red", "green", "blue")
g4 <- s.class(xy4, fac4, ppoints.cex = 1.5, chull = T, ellipseSize = 0, starSize = 0, ppolygon = list(border = 4:1, col = 1:4, lty = 1:4, lwd = 2, alpha = 0.4))

