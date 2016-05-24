library("symbolicDA")
data("cars",package="symbolicDA")
sdt<-cars
clust<-SClust(sdt, cl=3, iter=50)
print(clust)
