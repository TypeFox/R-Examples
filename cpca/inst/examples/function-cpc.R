require(plyr)
require(abind)

data(iris)

C <- daply(iris, "Species", function(x) cov(x[, -ncol(x)]))
C <- aperm(C, c(2, 3, 1)) # put the 1st dimension to the end

# default call
mod1 <- cpc(C)
round(mod1$CPC, 2)

# compute only first two CPCs
mod2 <- cpc(C, k = 2)
round(mod2$CPC, 2)

