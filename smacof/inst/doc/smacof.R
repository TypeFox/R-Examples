## ----echo=FALSE, hide=TRUE-----------------------------------------------
require(smacof)
require(RgoogleMaps)

## ------------------------------------------------------------------------
data(Guerry)
fit.guerry <- mds(Guerry)

## ----france-plot, eval=FALSE---------------------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fit.guerry)
#  theta <- 82*pi/180            ## degrees to radians
#  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
#  configs82 <- fit.guerry$conf %*% rot                ## rotated configurations
#  francemap1 <- GetMap(destfile="mypic1.png", zoom = 6, center = c(46.55, 3.05),
#                       maptype = "satellite")
#  PlotOnStaticMap(francemap1)
#  text(configs82*280, labels = rownames(configs82), col = "white", cex = 0.7)
#  par(op)

## ----france-plot1, echo=FALSE, fig.width=12, fig.height=6, dev='postscript'----
op <- par(mfrow = c(1,2))
plot(fit.guerry)
theta <- 82*pi/180            ## degrees to radians
rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
configs82 <- fit.guerry$conf %*% rot                ## rotated configurations
francemap1 <- GetMap(destfile="mypic1.png", zoom = 6, center = c(46.55, 3.05), 
                     maptype = "satellite")
PlotOnStaticMap(francemap1)
text(configs82*280, labels = rownames(configs82), col = "white", cex = 0.7)
par(op)

## ------------------------------------------------------------------------
ratings <- RockHard[,5:18]
rockdiss <- dist(t(ratings))
fit.rock <- mds(rockdiss, ndim = 3)
fit.rock

## ----judge-plot, eval=FALSE----------------------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fit.rock, plot.dim = c(1,2), main = "Configurations D1 vs. D2")
#  plot(fit.rock, plot.dim = c(1,3), main = "Configurations D1 vs. D3")
#  par(op)

## ----judge-plot1, echo=FALSE, fig.width=12, fig.height=6, dev='postscript'----
op <- par(mfrow = c(1,2))
plot(fit.rock, plot.dim = c(1,2), main = "Configurations D1 vs. D2")
plot(fit.rock, plot.dim = c(1,3), main = "Configurations D1 vs. D3")
par(op)

## ----shepard-plot, eval=FALSE--------------------------------------------
#  fit.interval <- mds(kinshipdelta, type = "interval")
#  fit.ordinal1 <- mds(kinshipdelta, type = "ordinal", ties = "primary")
#  fit.ordinal2 <- mds(kinshipdelta, type = "ordinal", ties = "secondary")
#  fit.spline <- mds(kinshipdelta, type = "mspline", spline.intKnots = 3,
#                    spline.degree = 2)
#  op <- par(mfrow = c(2,2))
#  plot(fit.interval, plot.type = "Shepard",
#       main = "Shepard Diagram (Interval MDS)", ylim = c(0.1, 1.7))
#  plot(fit.ordinal1, plot.type = "Shepard",
#       main = "Shepard Diagram (Ordinal MDS, Primary)", ylim = c(0.1, 1.7))
#  plot(fit.ordinal2, plot.type = "Shepard",
#       main = "Shepard Diagram (Ordinal MDS, Secondary)", ylim = c(0.1, 1.7))
#  plot(fit.spline, plot.type = "Shepard",
#       main = "Shepard Diagram (Spline MDS)", ylim = c(0.1, 1.7))
#  par(op)

## ----shepard-plot1, echo=FALSE, fig.width=8, fig.height=8, dev='postscript'----
fit.interval <- mds(kinshipdelta, type = "interval")
fit.ordinal1 <- mds(kinshipdelta, type = "ordinal", ties = "primary")
fit.ordinal2 <- mds(kinshipdelta, type = "ordinal", ties = "secondary")
fit.spline <- mds(kinshipdelta, type = "mspline", spline.intKnots = 3, 
                  spline.degree = 2)
op <- par(mfrow = c(2,2))
plot(fit.interval, plot.type = "Shepard", 
     main = "Shepard Diagram (Interval MDS)", ylim = c(0.1, 1.7))
plot(fit.ordinal1, plot.type = "Shepard", 
     main = "Shepard Diagram (Ordinal MDS, Primary)", ylim = c(0.1, 1.7))
plot(fit.ordinal2, plot.type = "Shepard", 
     main = "Shepard Diagram (Ordinal MDS, Secondary)", ylim = c(0.1, 1.7))
plot(fit.spline, plot.type = "Shepard", 
     main = "Shepard Diagram (Spline MDS)", ylim = c(0.1, 1.7))
par(op)

## ------------------------------------------------------------------------
LawlerD <- sim2diss(Lawler)
fitclas <- mds(LawlerD)
fitclas$stress
stressvec <- NULL
set.seed(123)
for(i in 1:20) {
  fitran <- mds(LawlerD, init = "random")
  stressvec[i] <- fitran$stress   
}
stressvec                          ## stress values

## ------------------------------------------------------------------------
stressvec <- randomstress(n = 9, ndim = 2, nrep = 500)
mean(stressvec)
fit <- mds(LawlerD)
fit$stress

## ----eval=FALSE----------------------------------------------------------
#  set.seed(1234)
#  res.perm <- permtest(fit, nrep = 1000, verbose = FALSE)
#  res.perm

## ----echo=FALSE----------------------------------------------------------
if(file.exists("resperm.rda")) load("resperm.rda") else {
set.seed(1234)
res.perm <- permtest(fit, nrep = 1000, verbose = FALSE)
}
res.perm

## ----perm-plot, eval=FALSE-----------------------------------------------
#  op <- par(mfrow = c(1,2))
#  hist(res.perm$stressvec, xlab = "Stress Values", main = "Histogram Permutations")
#  abline(v = quantile(res.perm$stressvec, c(0.025, 0.975)), col = "gray")
#  abline(v = fit$stress, col = "red", lwd = 2)
#  plot(res.perm)
#  par(op)

## ----perm-plot1, echo=FALSE, fig.height=5, fig.width=8, dev='postscript'----
op <- par(mfrow = c(1,2))
hist(res.perm$stressvec, xlab = "Stress Values", main = "Histogram Permutations")
abline(v = quantile(res.perm$stressvec, c(0.025, 0.975)), col = "gray")
abline(v = fit$stress, col = "red", lwd = 2)
plot(res.perm)
par(op)

## ------------------------------------------------------------------------
fit$spp

## ----spp-plot, eval=FALSE------------------------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fit, plot.type = "stressplot")
#  plot(fit, plot.type = "bubbleplot")
#  par(op)

## ----spp-plot1, echo=FALSE, fig.width = 9, fig.height = 5, dev='postscript'----
op <- par(mfrow = c(1,2))
plot(fit, plot.type = "stressplot")
plot(fit, plot.type = "bubbleplot")
par(op)

## ------------------------------------------------------------------------
ekmanD <- sim2diss(ekman, method = 1)
fit.basic <- mds(ekmanD, type = "ordinal")
fit.circ <- smacofSphere(ekmanD, type = "ordinal", verbose = FALSE)

## ----ekman-plot, eval=FALSE, echo=FALSE----------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fit.basic, main = "Unrestricted MDS")
#  plot(fit.circ, main = "Spherical MDS")
#  par(op)

## ----ekman-plot1, echo=FALSE, fig.width = 9, fig.height = 5, dev='postscript'----
op <- par(mfrow = c(1,2))
plot(fit.basic, main = "Unrestricted MDS")
plot(fit.circ, main = "Spherical MDS")
par(op)

## ------------------------------------------------------------------------
res.unc <- smacofSym(morse,  type = "ordinal")
res.parreg <- smacofConstraint(morse, type = "ordinal", ties = "primary", 
                               constraint = "linear", 
                               external = morsescales[,2:3], 
                               constraint.type = "ordinal", 
                               init = res.unc$conf)

## ----morse-plot, eval=FALSE----------------------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(res.unc, main = "Unconditional MDS")
#  plot(res.parreg, main = "Regional MDS")
#  par(op)

## ----morse-plot1, echo=FALSE, fig.width = 8, fig.height = 5, dev='postscript'----
op <- par(mfrow = c(1,2))
plot(res.unc, main = "Unconditional MDS")
plot(res.parreg, main = "Regional MDS")
par(op)

## ------------------------------------------------------------------------
ratings <- 11-RockHard[,5:18]            ## reverse ratings 
rownames(ratings) <- RockHard[,"Band"]
fit.rock <- unfolding(ratings)           ## 2D metric unfolding solution
fit.rock

## ----band-plot, eval=FALSE-----------------------------------------------
#  plot(fit.rock, label.conf.rows = list(label = FALSE))
#  best <- sort(rowMeans(ratings, na.rm = TRUE))[1:10]
#  worst <- sort(rowMeans(ratings, na.rm = TRUE), decreasing = TRUE)[1:10]
#  bestworst <- names(c(best, worst))
#  text(fit.rock$conf.row[bestworst,], labels = bestworst, cex = 0.8, pos = 3,
#       col = hcl(0, l = 50))

## ----band-plot1, echo=FALSE, fig.width = 7, fig.height = 7, dev='postscript'----
plot(fit.rock, label.conf.rows = list(label = FALSE))
best <- sort(rowMeans(ratings, na.rm = TRUE))[1:10]
worst <- sort(rowMeans(ratings, na.rm = TRUE), decreasing = TRUE)[1:10]
bestworst <- names(c(best, worst))
text(fit.rock$conf.row[bestworst,], labels = bestworst, cex = 0.8, pos = 3, 
     col = hcl(0, l = 50))

## ------------------------------------------------------------------------
PlatoD <- dist(t(Plato7))
fit.uni <- uniscale(PlatoD)
fit.uni

## ----d1-plot, eval=FALSE-------------------------------------------------
#  plot(fit.uni)

## ----d1-plot1, echo=FALSE, fig.width = 5, fig.height = 4, dev='postscript'----
plot(fit.uni)

## ------------------------------------------------------------------------
D <- as.matrix(kinshipdelta)[1:6, 1:6]
fit <- mds(D)                       ## MDS D --> conf
ifit <- inverseMDS(fit$conf)        ## inverse MDS conf --> D

## ----imds-plot, eval=FALSE-----------------------------------------------
#  op <- par(mfrow = c(3,3))
#  plot(fit, main = "Original MDS")
#  for (i in 1:length(ifit)) {
#    fit.i <- mds(ifit[[i]])
#    plot(fit.i, main = paste0("Inverse MDS (",i, ")"))
#  }
#  par(op)

## ----imds-plot1, echo=FALSE, fig.width = 8, fig.height = 8, dev='postscript'----
op <- par(mfrow = c(3,3))           
plot(fit, main = "Original MDS")
for (i in 1:length(ifit)) {
  fit.i <- mds(ifit[[i]])             
  plot(fit.i, main = paste0("Inverse MDS (",i, ")"))
}
par(op)

## ------------------------------------------------------------------------
eastD <- sim2diss(EW_eng$east)
attr(eastD, "Labels") <- abbreviate(attr(eastD, "Labels"))
fit.east <- mds(eastD, type = "ordinal")
westD <- sim2diss(EW_eng$west)
attr(westD, "Labels") <- abbreviate(attr(westD, "Labels"))
fit.west <- mds(westD, type = "ordinal", init = torgerson(eastD))

## ------------------------------------------------------------------------
fit.proc <- Procrustes(fit.east$conf, fit.west$conf)
fit.proc

## ----proc-plot, eval=FALSE-----------------------------------------------
#  op <- par(mfrow = c(2,2))
#  plot(fit.east, main = "MDS East Germany")
#  plot(fit.west, main = "MDS West Germany")
#  plot(fit.proc)
#  plot(fit.proc, plot.type = "transplot", length = 0.05)
#  par(op)

## ----proc-plot1, echo=FALSE, fig.width = 8, fig.height = 8, dev='postscript'----
op <- par(mfrow = c(2,2))
plot(fit.east, main = "MDS East Germany")
plot(fit.west, main = "MDS West Germany")
plot(fit.proc)
plot(fit.proc, plot.type = "transplot", length = 0.05)
par(op)

## ------------------------------------------------------------------------
fit.lawler <- mds(LawlerD, type = "interval")
jackfit <- jackknife(fit.lawler)
jackfit

## ----jack-plot, eval=FALSE-----------------------------------------------
#  plot(jackfit)

## ----jack-plot1, echo=FALSE, fig.width = 7, fig.height = 7, dev='postscript'----
plot(jackfit)

