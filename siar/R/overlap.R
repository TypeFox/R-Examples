overlap <- function (x1,y1,x2,y2,steps=5) {
# calculates the area of overlap between the shapes s1 and s2 which each
# have dimension dim(n,2) where each column are the x and y coordinates
# respectively.

E1 <- standard.ellipse(x1,y1,steps=steps)
E2 <- standard.ellipse(x2,y2,steps=steps)

ex1 <- E1$xSEAc[2:length(E1$xSEAc)]
ey1 <- E1$ySEAc[2:length(E1$xSEAc)]
ex2 <- E2$xSEAc[2:length(E2$xSEAc)]
ey2 <- E2$ySEAc[2:length(E2$xSEAc)]

e1 <- list(x=ex1,y=ey1)
e2 <- list(x=ex2,y=ey2)

out <- list()
out$overlap <- abs(overlap.xypolygon(e1,e2))
out$area1 <- abs(E1$SEAc)
out$area2 <- abs(E2$SEAc)

return(out)

}