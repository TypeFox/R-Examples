library("symbolicDA")
data("cars",package="symbolicDA")
r<- HINoV.SDA(cars, u=3, distance="U_2")
print(r$stopri)
plot(r$stopri[,2], xlab="Variable number", ylab="topri",
xaxt="n", type="b")
axis(1,at=c(1:max(r$stopri[,1])),labels=r$stopri[,1])
