## Exploring other representations of HE plots

iris.mod <- lm(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, data=iris)

vars <- c(3,1)
# get fitted values and residuals
iris.fit <- fitted(iris.mod)[,vars]
iris.res <- residuals(iris.mod)[,vars]

# translate residuals to grand mean
iris.ret <- scale(iris.res, center= -colMeans(iris[,vars]), scale=FALSE)

## ----------------------------------------------
## Explicitly showing residuals and fitted values
## ----------------------------------------------

heplot(iris.mod, var=vars,
xlab = "Petal Length in cm.", ylab = "Sepal length in cm.", size = "effect",
	main="HE plot: data ellipses of fitted values and residuals", 
	cex=1.25, cex.lab=1.25)

# take account of ranges in jittering fitted values
#apply(apply(iris.fit, 2, range), 2, diff)
#Petal.Length Sepal.Length 
#       4.090        1.582 

points(jitter(iris.fit[,1], amount=0, fac=2.5), jitter(iris.fit[,2], amount=0, fac=5), pch=16, col=rgb(0,0,1,.2))

legend(1.1, 6.8, c("fitted (jittered)", "residuals (centered)"), 
	col=c("blue", "red"), pch=16:15, cex=1.2)


## -----------------------------------------
## Contour plot of residuals & fitted values
## -----------------------------------------

if (require(KernSmooth)) {
  heplot(iris.mod, var=vars,
  xlab = "Petal Length in cm.", ylab = "Sepal length in cm.", size = "effect",
  	fill=c(TRUE, TRUE), fill.alpha=0.1,
  	cex=1.25, cex.lab=1.25)
  
  dens <- bkde2D(iris.ret, gridsize=c(61,61), bandwidth=0.3)
  contour(dens$x1, dens$x2, dens$fhat, levels=0.21, add=TRUE, col="brown", lwd=2)
  text(3.4, 5, "bkde2D(residuals)", cex=1.25, col="brown", pos=1)
  
  dens2 <- bkde2D(iris.fit, gridsize=c(61,61), bandwidth=0.3)
  contour(dens2$x1, dens2$x2, dens2$fhat, nlevels=2, add=TRUE, col="purple", lwd=2)
}


## --------------------
## bagplot of residuals
## --------------------

if (require(aplpack)) {
  heplot(iris.mod, var=vars,
  	xlab = "Petal Length in cm.", ylab = "Sepal length in cm.", size = "effect",
  	fill=c(TRUE, TRUE), fill.alpha=0.1, 
  #	level=0.5,
  	cex=1.25, cex.lab=1.25)
  
  iris.bag <- compute.bagplot(iris.ret[,1], iris.ret[,2])
  plot(iris.bag, add=TRUE, transparency=TRUE,
    show.outlier=FALSE, show.loophull=FALSE, col.baghull="#8B3333")  # brown3
}
