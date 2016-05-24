library(MSG)
data(ChinaLifeEdu)
x = ChinaLifeEdu
library(KernSmooth)
est = bkde2D(x, apply(x, 2, dpik))
contour(est$x1, est$x2, est$fhat, nlevels=15, col = 'darkgreen',
        vfont = c("sans serif", "plain"),
        xlab = "\u9884\u671F\u5BFF\u547D", ylab = "\u9AD8\u5B66\u5386\u4EBA\u6570")
points(x, pch = 20)
