### R code from vignette source 'plot3Drgl.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(plot3Drgl)
options(prompt = " ")
options(continue = "  ")
options(width=75)


###################################################
### code chunk number 2: plot3Drgl.Rnw:97-98
###################################################
persp3D(z = volcano, plot = FALSE)


###################################################
### code chunk number 3: plot3Drgl.Rnw:101-102
###################################################
args(plotrgl)


###################################################
### code chunk number 4: plot3Drgl.Rnw:108-109
###################################################
plotrgl(smooth = TRUE, lighting = TRUE)


###################################################
### code chunk number 5: volcano
###################################################
plotdev(shade = 0.1)


###################################################
### code chunk number 6: volcano
###################################################
plotdev(shade = 0.1)


###################################################
### code chunk number 7: plot3Drgl.Rnw:157-172
###################################################
attach(mtcars)
fit <- lm(mpg ~ wt + disp)

# predict values on regular xy grid
wt.pred <- seq(1.5, 5.5, length.out = 30)
disp.pred <- seq(71, 472, length.out = 30)
xy <- expand.grid(wt = wt.pred,
                  disp = disp.pred)

mpg.pred <- matrix (nrow = 30, ncol = 30,
    data = predict(fit, newdata = data.frame(xy),
    interval = "prediction"))

# fitted points for droplines to surface
fitpoints <- predict(fit)


###################################################
### code chunk number 8: fit
###################################################
scatter3D(z = mpg, x = wt, y = disp, colvar = abs(mpg - fitpoints),
      pch = 18, cex = 2, theta = 20, phi = 20, ticktype = "detailed",
      xlab = "wt", ylab = "disp", zlab = "mpg", main = "mtcars",
      clab = "error", zlim = c(5, 35),
      surf = list(x = wt.pred, y = disp.pred, z = mpg.pred,
                  facets = NA, border = "black", fit = fitpoints)
      )


###################################################
### code chunk number 9: fit
###################################################
scatter3D(z = mpg, x = wt, y = disp, colvar = abs(mpg - fitpoints),
      pch = 18, cex = 2, theta = 20, phi = 20, ticktype = "detailed",
      xlab = "wt", ylab = "disp", zlab = "mpg", main = "mtcars",
      clab = "error", zlim = c(5, 35),
      surf = list(x = wt.pred, y = disp.pred, z = mpg.pred,
                  facets = NA, border = "black", fit = fitpoints)
      )


###################################################
### code chunk number 10: plot3Drgl.Rnw:193-194
###################################################
detach(mtcars)


###################################################
### code chunk number 11: plot3Drgl.Rnw:197-198
###################################################
plotrgl(new = FALSE)


###################################################
### code chunk number 12: plot3Drgl.Rnw:209-218
###################################################
 x <- y <- z <- seq(-2, 2, length.out = 15)
 xyz <- mesh(x, y, z)
 F <- with(xyz, log(x^2 + y^2 + z^2 +
                10*(x^2 + y^2) * (y^2 + z^2) ^2))

# three levels, transparency added
 isosurf3D(x, y, z, F, level = seq(0, 4, by = 2),
   col = c("red", "blue", "yellow"),
   clab = "F", alpha = 0.2, plot = FALSE)


###################################################
### code chunk number 13: iso
###################################################
plotdev()


###################################################
### code chunk number 14: iso
###################################################
plotdev()


###################################################
### code chunk number 15: plot3Drgl.Rnw:232-233
###################################################
plotrgl(new = FALSE, lighting = TRUE)


