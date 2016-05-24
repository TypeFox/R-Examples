### R code from vignette source 'tclust.rnw'

###################################################
### code chunk number 1: chunk1
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

library (tclust)

                 ## standard margins
mmar <- c (5.1, 4.1, 4.1, 1.1)
og <- grey (0)   ## color for outliers


################
##  Figure 1  ##
################

set.seed (100)

data (geyser2)
clus <- tkmeans (geyser2, k = 3, alpha=0.03)
plot (clus, col = c (og, 2:4), tol.lwd = 1, tol.lty = 2)




###################################################
### code chunk number 2: chunk2
###################################################

################
##  Figure 2  ##
################

data (M5data)
cl <- M5data[, "cluster"]
plot (M5data[, 1:2], col = cl + 1, pch = cl + 1, main = "The M5data data set")



###################################################
### code chunk number 3: chunk3
###################################################

################
##  Figure 3  ##
################

data (M5data)
x <- M5data[, 1:2]
set.seed (100)

  ## applying tclust with restr.fact = 1, restr= "eigen"
res.a <- tclust (x, k = 3, alpha=0.1, restr.fact = 1, restr= "eigen",
                 equal.weights = TRUE, warnings = 1)

  ## applying tclust with restr= "sigma"
res.b <- tclust (x, k = 3, alpha=0.1, restr.fact = 1, restr= "sigma",
                 equal.weights = TRUE)

  ## applying tclust with restr.fact = 1, restr= "deter"
res.c <- tclust (x, k = 3, alpha=0.1, restr.fact = 1, restr= "deter",
                 equal.weights = TRUE, warnings = 1)

  ## applying tclust with restr.fact = 50, restr= "eigen"
res.d <- tclust (x, k = 3, alpha=0.1, restr.fact = 50, restr= "eigen",
                 equal.weights = FALSE)

op <- par (mfrow = c (2, 2), mar = c (2.1, 2.1, 3.6, 1.1))
col <- c (og, 4, 2, 3)
pch <- c (1, 4, 2, 3)

  ## plotting results..
plot (res.a, col = col, pch = pch, tol.lty = 2, xlab = "",
     ylab = "", main.pre ="(a)", main = "/r")
plot (res.b, col = col, pch = pch, tol.lty = 2, xlab = "",
     ylab = "", main.pre ="(b)", main = "/r")
plot (res.c, col = col, pch = pch, tol.lty = 2, xlab = "",
     ylab = "", main.pre ="(c)", main = "/r")
plot (res.d, col = col, pch = pch, tol.lty = 2, xlab = "",
     ylab = "", main.pre ="(d)", main = "/r")
par (op)




###################################################
### code chunk number 4: chunk4
###################################################

################
##  Figure 4  ##
################

  ## data generation
set.seed (10)
x <- rbind (rmvnorm (200, c (0, 0), diag (2)),
            rmvnorm (200, c (5, 0), diag (2)))

  ## applying tclust
clus.4 <- tclust (x, k = 3, alpha = 0.0, restr.fact = 1, warnings = 0)

  ## plotting results..
plot (clus.4, xlim = c (-3, 9), ylim = c (-6, 6))




###################################################
### code chunk number 5: chunk5
###################################################
################
##  Figure 5  ##
################

library (cluster)
data (geyser2)

set.seed (0)

op <- par (mfrow = c (1, 3), mar = c (3.1, 3.1, 3.1, 2.1))

xlab <- "" #names (geyser2)[1]
ylab <- "" #names (geyser2)[2]

  ## applying PAM on the geyser2 data set
clus <- pam (geyser2, 3)
plot (geyser2, xlab = xlab, ylab = ylab, col = clus$clustering + 1,
      main = "(a) PAM with Geyser")
points (clus$medoids, pch = "X", cex = 1.4)

  ## modifying the geyser2 data set
geyser2.mod <- geyser2
idx <- c (16,  21, 36, 171, 236, 265)
geyser2.mod[idx,] <- geyser2[idx,] - 12

  ## applying PAM on the modified geyser2 data set
clus <- pam (geyser2.mod, 3)
plot (geyser2.mod, xlab = xlab, ylab = ylab, col = clus$clustering + 1,
      main = "(b) PAM with modified Geyser")

  ## applying tkmeans on the modified geyser2 data set
clus <- tkmeans (geyser2.mod, k = 3, alpha = 0.03)
plot (clus, xlab = xlab, ylab = ylab,
      main = "(c) TCLUST with modified Geyser", sub = "")
par (op)



###################################################
### code chunk number 6: chunk6
###################################################
################
##  Figure 6  ##
################

fig.elt <- list ()

library (sn)

op <- par (mfrow = c (2, 3), mar = c (3.1, 2, 3.1, 2))

set.seed (3)

n <- 200
xlab <- "" #"x1"
ylab <- "" #"x2"

  ## data generation: Mixture of elliptical t's
xi <- c (0.5, -1)
Omega <- diag (2) / 2 + 0.5
alpha <- c (0, 0)
rnd1 <- rmst (n, xi, Omega, alpha, nu = 2)
xi <- c (6, -3)
Omega <- matrix (c (4, 1, 1, 2), nrow = 2, ncol = 2)
alpha <- c (0, 0)
rnd2 <- rmst (n, xi, Omega, alpha, nu = 3)
xi <- c (-3, 3)
Omega <- matrix (c (2, 1, 1, 4), nrow = 2, ncol = 2)
alpha <- c (0, 0)
rnd3 <- rmst (n, xi, Omega, alpha, nu = 4)
rnd.1 <- rbind (rnd1, rnd2, rnd3)
real.clus1 <- c (rep (1, n), rep (2, n), rep (3, n))

ylim <- c (-25, 15)
  ## Plotting the real clusters
plot (rnd.1, xlim = c (-15, 80), ylim = ylim, xlab = xlab, ylab = ylab,
      col = real.clus1 + 1, pch = real.clus1 + 1)
title ("(a) Elliptical t components", line = 1.6)

  ## applying tclust without trimming
clus <- tclust (rnd.1, k = 3, alpha = 0.0)
plot (clus, xlim = c (-15, 80), ylim = ylim, xlab = xlab, ylab = ylab,
      main.pre = "(b)")

 ## applying tclust with trimming (5%)
clus <- tclust (rnd.1, k = 3, alpha = 0.05)
plot (clus, xlim = c (-15, 80), ylim = ylim, xlab = xlab, ylab = ylab,
      main.pre = "(c)")

  ## data generation: Mixture of non-elliptical t's
xi <- c (0.5, -1)
Omega <- diag (2) / 2 + 0.5
alpha <- c (2, 80)
rnd1 <- rmst (n, xi, Omega, alpha, nu = 2)
xi <- c (6, -3)
Omega <- matrix (c (4, 1, 1, 2), nrow = 2, ncol = 2)
alpha <- c (2, 2)
rnd2 <- rmst (n, xi, Omega, alpha, nu = 3)
xi <- c (-3, 3)
Omega <- matrix (c (2, 1, 1, 4), nrow = 2, ncol = 2)
alpha <- c (4, 4)
rnd3 <- rmst (n, xi, Omega, alpha, nu = 4)

rnd.2 <- rbind (rnd1, rnd2, rnd3)
real.clus2 <- c (rep (1, n), rep (2, n), rep(3, n))

ylim <- c (-10, 20)
  ## Plotting the real clusters
plot (rnd.2, xlim = c (-10, 20), ylim = ylim, xlab = xlab, ylab = ylab,
      col = real.clus2 + 1, pch=real.clus2 + 1)
title ("(d) Non-elliptical t components", line = 1.6)

 ## applying tclust without trimming
clus <- tclust (rnd.2, k = 3, alpha = 0.0)
plot (clus, xlim = c (-10, 20), ylim = ylim, xlab = xlab, ylab = ylab,
      main.pre = "(e)")

 ## applying tclust with trimming (5%)
clus <- tclust (rnd.2, k = 3, alpha = 0.05)
plot (clus, xlim = c (-10, 20), ylim = ylim, xlab = xlab, ylab = ylab,
      main.pre = "(f)")

par (op)




###################################################
### code chunk number 7: chunk7
###################################################

################
##  Figure 7  ##
################

library (mclust)

  ## function for plotting mclust-objects in tclust-style.
mclustplottcluststyle2d <- function (x, clus, clus.perm, tol.lty = 3,
                                     size = sqrt (qchisq (0.95, 2)), ...)
{
  if (missing (clus.perm))
    cuse <- clus$classification
  else
    cuse <- clus.perm [clus$classification + 1]

  plot (x, col = cuse + 1, pch = cuse + 1, ...)

  k <- ncol (clus$parameters$mean)

  for (i in 1:k)
    tclust:::.doEllipses (cov = clus$parameters$variance$sigma[,, i],
                          center = clus$parameters$mean[, i], size = size,
                          lty = tol.lty)
}

op <- par (mfrow = c (2, 2), mar = c (3.1, 3.1, 3.1, 2.1))

  ## example 1 - uniformly distributed noise

  ## data generation

set.seed (0)
x <- rbind (
            rmvnorm (360, c (0.0,  0), matrix (c ( 1,  0,  0, 1), ncol = 2)),
            rmvnorm (540, c (5.0, 10), matrix (c ( 6, -2, -2, 6), ncol = 2)),
            cbind (runif (100, -5, 15), runif (100, -10, 20)))

  ## applying mclust
noiseInit <- sample (c (TRUE, FALSE), size = nrow (x),
                     replace = TRUE, prob = c (0.1, 0.9))
clus <- Mclust (x, initialization = list (noise = noiseInit), G = 2)
mclustplottcluststyle2d (x, clus, c (0, 2, 1), xlab = "", ylab = "",
                         main = "(a) mclust", tol.lty = 2)

  ##  applying tclust
clus <- tclust (x, k = 2, alpha = 0.1)
plot (clus, xlab = "", ylab = "", main = "(b) tclust", sub = "", tol.lty = 2)


  ## example 2 - structured noise

  ## data generation
set.seed (0)
v <- runif (100, -2 * pi, 2 * pi)
noise <- cbind (100 + 25 * sin (v), 10 + 5 * v)

x <- rbind (
            rmvnorm (360, c (0.0,  0), matrix (c (1,  0,  0, 1), ncol = 2)),
            rmvnorm (540, c (5.0, 10), matrix (c (6, -2, -2, 6), ncol = 2)),
            noise)

  ##  applying mclust
noiseInit <- sample (c (TRUE, FALSE), size = nrow (x),
                     replace = TRUE, prob = c (0.1, 0.9))
clus <- Mclust (x, initialization = list (noise = noiseInit), G = 2)
mclustplottcluststyle2d (x, clus, c (0, 2, 1), xlab = "", ylab = "",
                         main = "(c) mclust", tol.lty = 2)

  ##  applying tclust
clus <- tclust (x, k = 2, alpha = 0.1)
plot (clus, xlab = "", ylab = "", main = "(d) tclust", sub = "", tol.lty = 2)

par (op)




###################################################
### code chunk number 8: chunk8
###################################################

################
##  Figure 8  ##
################

  ## data generation
set.seed (100)
mixt <- rbind (rmvnorm (360, c (  0,  0), matrix (c (1,  0,  0,  1), ncol = 2)),
               rmvnorm (540, c (  5, 10), matrix (c (6, -2, -2,  6), ncol = 2)),
               rmvnorm (100, c (2.5,  5), matrix (c (50, 0,  0, 50), ncol = 2)))

  ## applying tclust
set.seed (100)
clus.1 <- tclust (mixt, k = 3, alpha =  0.0, restr.fact = 50, warnings = 2)
clus.2 <- tclust (mixt, k = 2, alpha = 0.05, restr.fact =  8, warnings = 2)

  ## plotting results
op <- par (mfrow = c (1, 2))
plot(clus.1, by.clust = TRUE, col = c (og, 2, 3, 4), pch = c (1, 2, 3, 4),
     tol.lty = 2, main.pre = "(a)")
plot(clus.2, by.clust = TRUE, col = c (og, 2, 3), pch = c (1, 2, 3),
     tol.lty = 2, main.pre = "(b)")
par (op)




###################################################
### code chunk number 9: chunk9
###################################################

################
##  Figure 9  ##
################

  ## computing ctlcurves
ctl <- ctlcurves (mixt, restr.fact = 50, alpha = seq (0, 0.2, by = 0.05))

op <- par (mfrow = c (1, 1))
plot (ctl)
par (op)



###################################################
### code chunk number 10: chunk10
###################################################
#################
##  Figure 10  ##
#################

data (geyser2)
op <- par (mfrow = c (1, 2))
set.seed (10)

  ## creating a one dimensional data "matrix"
geyser1 <- geyser2[, 1, drop = FALSE]
  ## applying tkmeans
plot (tkmeans (geyser1, k = 2, alpha = 0.03), jitter = TRUE, tol.lwd = 2,
      main.pre = "(a)")

  ## adding a random dimension to geyser2
geyser3 <- cbind (geyser2, rnorm (nrow (geyser2)))
  ## applying tkmeans
plot (tkmeans (geyser3, k = 3, alpha = 0.03), main.pre = "(b)")

par (op)



###################################################
### code chunk number 11: chunk11
###################################################

#################
##  Figure 11  ##
#################

  ## data generation
set.seed (100)
mixt2 <- rbind (
  rmvnorm (360, c (0,   0), matrix (c (1,  0,  0,  1), ncol = 2)),
  rmvnorm (540, c (5,  10), matrix (c (6, -2, -2,  6), ncol = 2)),
  rmvnorm (100, c (2.5, 5), matrix (c (50, 0,  0, 50), ncol = 2))
)

  ## applying tclust
set.seed (100)
clus.w <- tclust (mixt2, k = 3, alpha = 0.1, restr.fact = 1,
                  equal.weights = TRUE, warnings = 1)

  ## applying DiscrFact
discr.clus.w <- DiscrFact (clus.w)

  ## plotting results
op <- par (mfrow = c (1, 3), mar = mmar)

plot (clus.w, col = c (og, 2, 3, 4), tol.lty = 2, main.pre = "(a)")
plot_DiscrFact_p2 (discr.clus.w, xlim = c (-70, 0), main.pre = "(b)")
plot_DiscrFact_p3 (discr.clus.w, tol.lty = 2, main.pre = "(c)")

par (op)


###################################################
### code chunk number 12: chunk12
###################################################
#################
##  Figure 12  ##
#################

data (swissbank)
set.seed (0)
fig6.ctl <- ctlcurves (swissbank, k = 1:4, alpha = seq (0, 0.3, by = 0.025),
                       restr.fact = 50, iter.max = 100, nstart = 100)

op <- par (mfrow = c (1, 1), mar = mmar)
plot (fig6.ctl)

par (op)


###################################################
### code chunk number 13: chunk13
###################################################
#################
##  Figure 13  ##
#################

fig7.clus <- tclust (swissbank, k = 2, alpha = .1, restr.fact = 50)
fig7.discrfact <- DiscrFact (fig7.clus, threshold = .000125)

op <- par (mfrow = c (1, 3), mar = mmar)
plot (fig7.discrfact, enum = TRUE)
par (op)


###################################################
### code chunk number 14: chunk14
###################################################

#################
##  Figure 14  ##
#################


clus <- tclust (swissbank, k = 2, alpha=.1, restr.fact=50)
discrfact <- DiscrFact (clus)
pch <- c (rep ("G", 100), rep ("F", 100))
condition <- discrfact$assignfact > log (.000125)
cl <- clus$cluster
data (swissbank)
op <- par (mfrow = c (1, 3), mar = mmar)

xlab <- "Distance of the inner frame to lower border"
ylab <- "Length of the diagonal"

plot (swissbank[, 4], swissbank[, 6], col = "darkgrey", pch = pch,
      main = "(a) Cluster1", xlab = xlab, ylab = ylab)

cl1 <- cl == 1
points (swissbank[cl1, 4], swissbank[cl1, 6], pch = pch[cl1], col = 2)
idx <- (cl1) & condition
points (swissbank[idx, 4], swissbank[idx, 6], pch = 1, cex = 4, col = "blue")

plot (swissbank[, 4], swissbank[, 6], col = "darkgrey", pch = pch,
      main = "(b) Cluster2", xlab = xlab, ylab = ylab)
cl2 <- cl == 2
points (swissbank[cl2, 4], swissbank[cl2, 6], pch = pch[cl2], col = 3)
idx <- (cl2) & condition
points (swissbank[idx, 4], swissbank[idx, 6], pch = 1, cex = 4, col = "blue")

cl0 <- cl == 0
plot (swissbank[, 4], swissbank[, 6], col = "darkgrey", pch = pch,
      main = "(c) Trimmed", xlab = xlab, ylab = ylab)

points (swissbank[cl0, 4], swissbank[cl0, 6], pch = pch[cl0])
idx <- (cl0) & condition
points (swissbank[idx, 4], swissbank[idx, 6], pch = 1, cex = 4, col = "blue")
par (op)




