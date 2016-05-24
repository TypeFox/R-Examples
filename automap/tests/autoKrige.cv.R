library(automap)

data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

kr.cv = autoKrige.cv(log(zinc)~1, meuse, model = c("Exp"))
kr_dist.cv = autoKrige.cv(log(zinc)~sqrt(dist), meuse, 
       model = c("Exp"))
kr_dist_ffreq.cv = autoKrige.cv(log(zinc)~sqrt(dist)+ffreq, 
       meuse, model = c("Exp"))

summary(kr.cv)
summary(kr_dist.cv)
summary(kr_dist_ffreq.cv)

compare.cv(kr.cv, kr_dist.cv, kr_dist_ffreq.cv)