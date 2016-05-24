## This file checks if the analytical gradient matches the numerical
## one for the max-stable models
library(SpatialExtremes)
n.site <- 20
n.obs <- 50

set.seed(1)
coord <- matrix(runif(2 * n.site, 0, 10), n.site)
data <- rmaxstab(n.obs, coord, "gauss", cov11 = 3, cov12 = 0, cov22 = 3)

M1 <- fitmaxstab(log(data), coord, "gauss", y ~ 1, y ~ 1, y ~ 1,
                 check.grad = TRUE, iterlim = 2)

M2 <- fitmaxstab(log(data), coord, "whitmat", y ~ 1, y ~ 1, y ~ 1,
                 check.grad = TRUE, iterlim = 2)

M3 <- fitmaxstab(log(data), coord, "gpowexp", y ~ 1, y ~ 1, y ~ 1,
                 check.grad = TRUE, iterlim = 2)

M4 <- fitmaxstab(log(data), coord, "brown", y ~ 1, y ~ 1, y ~ 1,
                 check.grad = TRUE, iterlim = 2)
