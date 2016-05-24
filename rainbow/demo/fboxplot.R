## This demo considers functional data visualization, including functional 
## bagplot and functional HDR boxplot

library(rainbow)

fboxplot(data = ElNino, plot.type = "functional", type = "bag")
fboxplot(data = ElNino, plot.type = "bivariate", type = "bag")
fboxplot(data = ElNino, plot.type = "functional", type = "hdr", alpha = c(0.07, 0.5))
fboxplot(data = ElNino, plot.type = "bivariate", type = "hdr", alpha = c(0.07, 0.5))

