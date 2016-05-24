library("symbolicDA")
data("cars",package="symbolicDA")
dist<-dist.SDA(cars, type="U_3", gamma=0.3, power=2)
print(dist)
