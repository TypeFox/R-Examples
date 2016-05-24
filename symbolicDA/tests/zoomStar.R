library(symbolicDA) 
# Example 1
data("cars",package="symbolicDA")
sdt<-cars
zoomStar(sdt, j=12)

# Example 2
data("cars",package="symbolicDA")
sdt<-cars
variables<-as.matrix(sdt$variables)
indivN<-as.matrix(sdt$indivN)
dist<-as.matrix(dist.SDA(sdt))
classes<-DClust(dist, cl=5, iter=100)
for(i in 1:max(classes)){
  getOption("device")()  
  zoomStar(sdt, .medoid(dist, classes, i))}
