source("../R/distance3d.r")

detectGlobal(dist3d)
x <- matrix(c(1,2,3),1)
y <- matrix(c(0,0,0),1)
dist3d(x,y)
dist3d(x,y,method="maximum")


x <- matrix(c(1,2,3,2,3,4),2,byrow=TRUE)
y <- matrix(c(0,0,0,3,3,3),2,byrow=TRUE)
dist3d(x,y,method="manhattan")
dist3d(x,y)
dist3d(x,y,method="minkowski",p=3)
dist3d(x,y,method="maximum")

