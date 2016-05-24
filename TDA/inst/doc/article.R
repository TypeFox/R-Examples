### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: initial_settings
###################################################
options(stringsAsFactors = FALSE)
options(width = 65)
options(prompt = "> ")


###################################################
### code chunk number 2: eq1
###################################################
library("TDA")
X <- circleUnif(400)

Xlim <- c(-1.6, 1.6);  Ylim <- c(-1.7, 1.7);  by <- 0.065

Xseq <- seq(Xlim[1], Xlim[2], by = by)
Yseq <- seq(Ylim[1], Ylim[2], by = by)
Grid <- expand.grid(Xseq, Yseq)


###################################################
### code chunk number 3: eq2
###################################################
distance <- distFct(X = X, Grid = Grid)


###################################################
### code chunk number 4: eq3
###################################################
m0 <- 0.1
DTM <- dtm(X = X, Grid = Grid, m0 = m0)


###################################################
### code chunk number 5: eq4
###################################################
k <- 60
kNN <- knnDE(X = X, Grid = Grid, k = k)


###################################################
### code chunk number 6: eq5
###################################################
h <- 0.3
KDE <- kde(X = X, Grid = Grid, h = h)


###################################################
### code chunk number 7: eq6
###################################################
h <- 0.3
Kdist <- kernelDist(X = X, Grid = Grid, h = h)


###################################################
### code chunk number 8: eq7 (eval = FALSE)
###################################################
## persp(Xseq, Yseq,
##       matrix(KDE, ncol = length(Yseq), nrow = length(Xseq)), xlab = "",
##       ylab = "", zlab = "", theta = -20, phi = 35, ltheta = 50,
##       col = 2, border = NA, main = "KDE", d = 0.5, scale = FALSE,
##       expand = 3, shade = 0.9)


###################################################
### code chunk number 9: eq8
###################################################
par(mfrow = c(2, 3))
par(mai = c(0.6, 0.25, 0.3, 0.1))
plot(X, pch = 16, cex = 0.6, xlab = "", ylab = "", main = "Sample X")

par(mai = c(0.1, 0.15, 0.1, 0.1))
persp(Xseq, Yseq,
      matrix(distance, ncol = length(Yseq), nrow = length(Xseq)),
      xlab = "", ylab = "", zlab = "", theta = -20, phi = 35,
      ltheta = 50, col = 2, border = NA, main = "Distance Function",
      d = 0.5, scale = FALSE, expand = 1, shade = 0.9)
persp(Xseq, Yseq,
      matrix(DTM, ncol = length(Yseq), nrow = length(Xseq)), xlab = "",
      ylab = "", zlab = "", theta = -20, phi = 35, ltheta = 50,
      col = 2, border = NA, main = "DTM", d = 0.5, scale = FALSE,
      expand = 1, shade = 0.9)
persp(Xseq, Yseq,
      matrix(kNN, ncol = length(Yseq), nrow = length(Xseq)), xlab = "",
      ylab = "", zlab = "", theta = -20, phi = 35, ltheta = 50,
      col = 2, border = NA, main = "kNN", d = 0.5, scale = FALSE,
      expand = 3, shade = 0.9)
persp(Xseq, Yseq,
      matrix(KDE, ncol = length(Yseq), nrow = length(Xseq)), xlab = "",
      ylab = "", zlab = "", theta = -20, phi = 35, ltheta = 50,
      col = 2, border = NA, main = "KDE", d = 0.5, scale = FALSE,
      expand = 3, shade = 0.9)
persp(Xseq, Yseq,
      matrix(Kdist, ncol = length(Yseq), nrow = length(Xseq)),
      xlab = "", ylab = "", zlab = "", theta = -20, phi = 35,
      ltheta = 50, col = 2, border = NA, main = "Kernel Distance",
      d = 0.5, scale = FALSE, expand = 5, shade = 0.9)


###################################################
### code chunk number 10: eq9
###################################################
band <- bootstrapBand(X = X, FUN = kde, Grid = Grid, B = 100,
                      parallel = FALSE, alpha = 0.1, h = h)


###################################################
### code chunk number 11: eq10
###################################################
posYgrid <- which(Grid[, 2] > 0)
posYseq <- which(Yseq > 0)

par(mfrow = c(1, 2))
par(mai = c(0.5, 0.25, 0.3, 0.1))

persp(Xseq, Yseq, matrix(band[["band"]][, 1], ncol = length(Yseq),
      nrow = length(Xseq)),
      zlim = c(0, max(band[["band"]][posYgrid, 2])),
      ylim = range(Yseq), xlab = "", ylab = "", zlab="",theta = 0,
      phi = 25, ltheta = 50, col = "pink", border = NA, d = 0.5,
      scale = FALSE, expand = 3, shade = 0.9, box = FALSE)


persp(Xseq, Yseq[posYseq], matrix(band[["band"]][posYgrid, 1],
      ncol = length(Yseq[posYseq]), nrow = length(Xseq)),
      zlim = c(0, max(band[["band"]][posYgrid, 2])),
      ylim = c(-0.5, 1.7), xlab = "", ylab = "", zlab = "", theta = 0,
      phi = 25, ltheta = 50, col = "pink", border = NA, d = 0.5,
      scale = FALSE, expand = 3, shade = 0.9, box = FALSE)

par(new = TRUE)

persp(Xseq, Yseq[posYseq], matrix(band[["fun"]][posYgrid],
      ncol = length(Yseq[posYseq]), nrow = length(Xseq)),
      zlim = c(0, max(band[["band"]][posYgrid, 2])),
      ylim = c(-0.5, 1.7), xlab = "", ylab = "", zlab = "", theta = 0,
      phi = 25, ltheta = 50, col = "red", border = NA, d = 0.5,
      scale = FALSE, expand = 3, shade = 0.9, box = FALSE)

par(new = TRUE)

persp(Xseq, Yseq[posYseq], matrix(band[["band"]][posYgrid, 2],
      ncol = length(Yseq[posYseq]), nrow = length(Xseq)),
      zlim = c(0, max(band[["band"]][posYgrid, 2])),
      ylim = c(-0.5, 1.7), xlab = "", ylab = "", zlab = "", theta = 0,
      phi = 25, ltheta = 50, col = "pink", border = NA, main = "",
      d = 0.5, scale = FALSE, expand = 3, shade = 0.9, box = FALSE)


###################################################
### code chunk number 12: eq11
###################################################
Diag <- gridDiag(X = X, FUN = kde, h = 0.3, lim = cbind(Xlim, Ylim),
                 by = by, sublevel = FALSE, library = "Dionysus",
                 printProgress = FALSE)


###################################################
### code chunk number 13: eq11b (eval = FALSE)
###################################################
## plot(Diag[["diagram"]], band = 2 * band[["width"]],
##      main = "KDE Diagram")


###################################################
### code chunk number 14: eq11b
###################################################
par(mfrow = c(1, 3))
par(mai = c(0.3, 0.25, 0.3, 0.1))
plot(X, pch = 16, cex = 0.6, xlab = "", ylab = "", main = "Sample X")
persp(Xseq, Yseq,
      matrix(KDE, ncol = length(Yseq), nrow = length(Xseq)),
      xlab = "", ylab = "", zlab = "", theta = -20, phi = 35,
      ltheta = 50, col = 2, border = NA, main = "KDE", d = 0.5,
      scale = FALSE, expand = 3, shade = 0.9)
par(mai = c(0.55, 0.6, 0.3, 0.1))
plot(Diag[["diagram"]], band = 2 * band[["width"]],
     main = "KDE Diagram")


###################################################
### code chunk number 15: eq11c
###################################################
par(mfrow = c(1, 2), mai = c(0.8, 0.8, 0.3, 0.1))
plot(Diag[["diagram"]], rotated = TRUE, band = band[["width"]],
     main = "Rotated Diagram")
plot(Diag[["diagram"]], barcode = TRUE,  main = "Barcode")


###################################################
### code chunk number 16: eq12
###################################################
Circle1 <- circleUnif(60)
Circle2 <- circleUnif(60, r = 2) + 3
Circles <- rbind(Circle1, Circle2)


###################################################
### code chunk number 17: eq12b
###################################################
maxscale <- 5        # limit of the filtration
maxdimension <- 1    # components and loops


###################################################
### code chunk number 18: eq12c
###################################################
Diag <- ripsDiag(X = Circles, maxdimension, maxscale,
                 library = "GUDHI", printProgress = FALSE)


###################################################
### code chunk number 19: eq12d
###################################################
par(mfrow = c(1, 2), mai=c(0.8, 0.8, 0.3, 0.3))
plot(Circles, pch = 16, xlab = "",ylab = "")
plot(Diag[["diagram"]])


###################################################
### code chunk number 20: eq13
###################################################
Diag1 <- ripsDiag(Circle1, maxdimension = 1, maxscale = 5)
Diag2 <- ripsDiag(Circle2, maxdimension = 1, maxscale = 5)


###################################################
### code chunk number 21: eq13b
###################################################
print(bottleneck(Diag1[["diagram"]], Diag2[["diagram"]],
                 dimension = 1))
print(wasserstein(Diag1[["diagram"]], Diag2[["diagram"]], p = 2,
                  dimension = 1))


###################################################
### code chunk number 22: eq14a
###################################################
PlotTriangles <- function(left, right) {
  n <- length(left)
  X <- (left + right) / 2
  Y <- (right - left) / 2
  xlimit <- c(0, max(X + Y) * 1.2)
  ylimit <- c(0, max(Y) * 1.2)
  plot(X, Y, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       xlim = xlimit, ylim = ylimit, frame.plot = FALSE,
       main = "Triangles")
  axis(1)
  axis(2)
	mtext("(Birth + Death) / 2", side = 1, line = 2.5, cex = 1)
	mtext("(Death - Birth) / 2", side = 2, line = 2.5, cex = 1)
	#polygon(c(0, 0, ylimit[2]), c(0, ylimit[2], ylimit[2]),
  #        col = "gray87", border = NA)
  for (i in seq_len(n)) {
    a <- c(0, left[i], (left[i] + right[i]) / 2, right[i], 0)
    b <- c(0, 0, (right[i] - left[i]) / 2, 0, 0)
    lines(a, b, lwd = 2)
  }
  B <- max(X + Y)
  grid <- seq(0, B, length = 1000)
  maxfun <- rep(0, length(grid))
  tmp <- matrix(0, n, length(grid))
  for(i in seq_len(n)) {
    tmp[i, ] <- pmax(pmin(grid - left[i], right[i] - grid), 0)
  }
  land <- apply(tmp, 2, max)
  #lines(grid, land, lwd = 2, col = 4)
	points(X, Y, type = "p", cex = 1, pch = 19, col = "pink")
  points(X, Y, type = "p")
}

par(mfrow = c(1, 3), mai = c(0.7, 0.7, 0.3, 0.3))

left <- c(0, 2, 2, 3.5, 5)
right <- c(2, 6, 3, 5, 8)

PlotTriangles(left, right)    

Diag1 <- cbind(rep(0, length(left)), left, right)
tseq <- seq(min(Diag1[, 2:3]), max(Diag1[, 2:3]), length = 500)
Land1 <- landscape(Diag1, dimension = 0, KK = 1, tseq)
Sil1 <- silhouette(Diag1, p = 1, dimension = 0, tseq)

X <- (left + right) / 2
Y <- (right - left) / 2
xlimit <- c(0, max(X + Y) * 1.2)
ylimit <- c(0, max(Y) * 1.2)

plot(tseq, Land1, col = 4, type = "l", lwd = 2, xlim = xlimit,
     ylim = ylimit, axes = FALSE, main = "1st Landscape", xlab = "",
     ylab = "")
axis(1)
axis(2)

plot(tseq, Sil1, col = 4, type = "l", lwd = 2, xlim = xlimit,
     ylim = ylimit, axes = FALSE, main = "Silhouette p = 1", xlab = "",
     ylab = "")
axis(1)
axis(2)


###################################################
### code chunk number 23: eq14
###################################################
tseq <- seq(0, maxscale, length = 1000)   #domain
Land <- landscape(Diag[["diagram"]], dimension = 1, KK = 1, tseq)
Sil <- silhouette(Diag[["diagram"]], p = 1, dimension = 1, tseq)


###################################################
### code chunk number 24: eq14b
###################################################
par(mfrow = c(1, 2), mai = c(0.5, 0.45, 0.3, 0.3))
plot(tseq, Land, type = "l", lwd = 3, ylab = "",
     main = "1st Landscape, dim = 1", asp = 1, col = 2)
plot(tseq, Sil, type = "l", lwd = 3, ylab = "",
     main = "Silhouette(p = 1), dim = 1", asp = 1, col = 2)


###################################################
### code chunk number 25: eq15
###################################################
N <- 4000
XX1 <- circleUnif(N / 2)
XX2 <- circleUnif(N / 2, r = 2) + 3
X <- rbind(XX1, XX2)


###################################################
### code chunk number 26: eq15b
###################################################
m <- 80     # subsample size
n <- 10     # we will compute n landscapes using subsamples of size m
tseq <- seq(0, maxscale, length = 500)          #domain of landscapes

#here we store n Rips diags
Diags <- list()
#here we store n landscapes
Lands <- matrix(0, nrow = n, ncol = length(tseq))


###################################################
### code chunk number 27: eq15c
###################################################
for (i in seq_len(n)) {
  subX <- X[sample(seq_len(N), m), ]
  Diags[[i]] <- ripsDiag(subX, maxdimension = 1, maxscale = 5)
  Lands[i, ] <- landscape(Diags[[i]][["diagram"]], dimension = 1,
                          KK = 1, tseq)
}


###################################################
### code chunk number 28: eq15d
###################################################
bootLand <- multipBootstrap(Lands, B = 100, alpha = 0.05,
                            parallel = FALSE)


###################################################
### code chunk number 29: eq15e (eval = FALSE)
###################################################
## plot(tseq, bootLand[["mean"]], main = "Mean Landscape with 95% band")
## polygon(c(tseq, rev(tseq)),
##         c(bootLand[["band"]][, 1], rev(bootLand[["band"]][, 2])),
##         col = "pink")
## lines(tseq, bootLand[["mean"]], lwd = 2, col = 2)


###################################################
### code chunk number 30: eq15f
###################################################
par(mfrow = c(1, 2))
par(mai = c(0.5, 0.45, 0.3, 0.3))
plot(X, pch = 16, cex = 0.5, main = "Large Sample from Circles")
plot(tseq, bootLand[["mean"]], type = "l", lwd = 2, xlab = "",
     ylab = "", main = "Mean Landscape with 95% band",
     ylim = c(0, 1.2))
polygon(c(tseq, rev(tseq)),
        c(bootLand[["band"]][, 1], rev(bootLand[["band"]][, 2])),
        col = "pink")
lines(tseq, bootLand[["mean"]], lwd = 2, col = 2)


###################################################
### code chunk number 31: eq16
###################################################
XX1 <- circleUnif(600)
XX2 <- circleUnif(1000, r = 1.5) + 2.5
noise <- cbind(runif(80, -2, 5), runif(80, -2, 5))
X <- rbind(XX1, XX2, noise)

# Grid limits
Xlim <- c(-2, 5)
Ylim <- c(-2, 5)
by <- 0.2


###################################################
### code chunk number 32: eq16b
###################################################
parametersKDE <- seq(0.1, 0.6, by = 0.05)

B <- 50       # number of bootstrap iterations. Should be large.
alpha <- 0.1  # level of the confidence bands


###################################################
### code chunk number 33: eq16c
###################################################
maxKDE <- maxPersistence(kde, parametersKDE, X,
              lim = cbind(Xlim, Ylim), by = by, sublevel = FALSE,
              B = B, alpha = alpha, parallel = TRUE,
              printProgress = TRUE, bandFUN = "bootstrapBand")


###################################################
### code chunk number 34: eq16d
###################################################
print(summary(maxKDE))


###################################################
### code chunk number 35: eq16e
###################################################
par(mfrow = c(1, 2), mai = c(0.8, 0.8, 0.35, 0.3))
plot(X, pch = 16, cex = 0.5, main = "Two Circles")
plot(maxKDE, main = "Max Persistence - KDE")


###################################################
### code chunk number 36: eq18
###################################################
X1 <- cbind(rnorm(300, 1, .8), rnorm(300, 5, 0.8))
X2 <- cbind(rnorm(300, 3.5, .8), rnorm(300, 5, 0.8))
X3 <- cbind(rnorm(300, 6, 1), rnorm(300, 1, 1))
XX <- rbind(X1, X2, X3)


###################################################
### code chunk number 37: eq18b
###################################################
Tree <- clusterTree(XX, k = 100, density = "knn",
                    printProgress = FALSE)
TreeKDE <- clusterTree(XX, k = 100, h = 0.3, density = "kde",
                       printProgress = FALSE)


###################################################
### code chunk number 38: eq18c (eval = FALSE)
###################################################
## plot(Tree, type = "lambda", main = "lambda Tree (knn)")
## plot(Tree, type = "kappa", main = "kappa Tree (knn)")
## plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")
## plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")


###################################################
### code chunk number 39: eq18d
###################################################

par(mfrow = c(2,3))
par(mai = c(0.25,0.35,0.3,0.3))

plot(XX, pch = 16, cex = 0.6, main = "Data")

# plot lambda trees
plot(Tree, type = "lambda", main = "lambda Tree (knn)")
plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")

# plot clusters
plot(XX, pch = 19, cex = 0.6, main = "cluster labels (knn)")
for (i in Tree[["id"]]) {
	points(matrix(XX[Tree[["DataPoints"]][[i]], ], ncol = 2), col = i,
         pch = 19, cex = 0.6)
}

#plot kappa trees
plot(Tree, type = "kappa", main = "kappa Tree (knn)")
plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")


