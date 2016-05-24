library("symbolicDA")
data("cars",package="symbolicDA")
y<-cars
cl<-SClust(y, 4, iter=150)
print(cl)
o<-cluster.Description.SDA(y, cl)
print(o)
