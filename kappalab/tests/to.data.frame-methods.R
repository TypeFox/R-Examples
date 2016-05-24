library(kappalab)

a <- Mobius.set.func(runif(16),4,4)
d <- to.data.frame(a)
write.table(d,"my.set.func.csv",sep="\t")
d2 <- read.delim("my.set.func.csv")
a2 <- Mobius.set.func(d2[[1]][3:18],4,4)
stopifnot(abs(a@data - a2@data) < 1e-6)



