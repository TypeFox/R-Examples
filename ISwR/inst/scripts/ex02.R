x <- 1:10
z <- append(x, 1.23, after=7) 
z
z <- c(x[1:7],1.23,x[8:10])
z
v <- 1.23; k <- 7
i <- seq(along=x)
z <- c(x[i <= k], v, x[i > k])
z
write.table(thuesen, file="foo.txt")
# edit the file
read.table("foo.txt", na.strings=".")
write.table(thuesen, file="foo.txt", na=".")
read.table("foo.txt", na.strings=".")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
