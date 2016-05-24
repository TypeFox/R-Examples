##test TP variant of svcm on 3d brain data; global grid search
library(svcm)
data(brain3d)
X <- matrix(c(0.5, 0.5,   0,   0, 0.5, 0.5,
                0,   0, 0.5, 0.5, 0.5, 0.5,
              0.5, 0.5, 0.5, 0.5,   0,   0,
                0,   0,   0,   0,   1,  -1,
                1,  -1,   0,   0,   0,   0,
                0,   0,   1,  -1,   0,   0), nrow = 6)
M3d <- svcm(brain3d, X, knots=c(3, 8, 4), deg=c(1, 1, 1),
          vsize=c(1.875, 1.875, 4.0), search=TRUE, type="TP",
          lambda.init=0.1, lower=-10, upper=0, ngrid=3)
