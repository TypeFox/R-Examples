# ex1.12.R
require(sde)
set.seed(123)
plot(BM())
plot(GBM(1,1,sqrt(0.5)))
plot(BBridge(0,-1))
