### Weighted and robust archetypal analysis: artificial toy data set
###
### Analysis used in 'Weighted and Robust Archetypal Analysis' by
### Manuel J. A. Eugster and Friedrich Leisch.

library('archetypes')



### Data set:

data('toy')

plot(toy, col = 1, pch = 19)



### Original archetypes:

set.seed(1986)
a <- archetypes(toy, 3)

a
parameters(a)
barplot(a, toy, percentiles = TRUE)


## Approximation of the convex hull:
xyplot(a, toy, chull = chull(toy))


## Approximtion of the data:
xyplot(a, toy, adata.show = TRUE)


## Residuals and RSS:
local({
  R <- residuals(a)

  opar <- par(mfrow = c(2, 1))
  plot(R[, 1], ylab = 'x residuals')
  plot(R[, 2], ylab = 'y residuals')
  par(opar)
})

plot(rss(a, type = 'single'), ylab = 'RSS (single)')


## Visualization of the iterations:
opar <- par(ask = TRUE)
movieplot(a, toy)
par(opar)



### Weighted archetypes:

w <- rep(1, nrow(toy))
w[toy[, 1] < 5] <- 0.3

set.seed(1234)
wa <- weightedArchetypes(toy, 3, weights = w)

wa
barplot(wa, toy, percentiles = TRUE)


## Data:
plot(toy, pch = 19, col = gray(1 - w))


## Weighted approximation of the convex hull:
xyplot(wa, toy)


## Weighted approximation of the data:
xyplot(wa, toy, adata.show = TRUE)


## Visualization of the iterations:
opar <- par(ask = TRUE)
movieplot(wa, toy)
par(opar)



### Data set with outliers:

toy.outlier <- function(n, mean, sigma) {
  require(mvtnorm)
  data(toy, package = 'archetypes')

  for ( i in seq(length = length(n)) )
    toy <- rbind(toy, rmvnorm(n[i], mean[[i]], sigma[[i]]))

  toy
}

set.seed(1234)
toy.o1 <- toy.outlier(5, mean = list(c(3, 28)),
                      sigma = list(diag(0.5, 2)))


plot(toy.o1, col = 1, pch = 19)



### Original archetypes:

set.seed(1234)
a1 <- archetypes(toy.o1, 3)

a1
xyplot(a1, toy.o1)



### Robust archetypes:

set.seed(1234)
ra <- robustArchetypes(toy.o1, 3)

ra
parameters(ra)
barplot(ra, toy.o1, percentiles = TRUE)


## Robust approximation of the convex hull:
xyplot(ra, toy.o1)


## Robust approximation of the data:
xyplot(ra, toy.o1, adata.show = TRUE)


## Panorama and RSS:
panorama(ra, toy.o1)
plot(rss(ra, type = 'single'), ylab = 'RSS (single)')


## Visualization of the iterations:
opar <- par(ask = TRUE)
movieplot(ra, toy.o1)
par(opar)
