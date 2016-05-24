library(wskm)

csv <- "ewkm01.csv"

ds <- read.csv("ewkm01.csv")

set.seed(42)
km <- ewkm(ds, 10)
print(km)
cat(km$iterations, km$restarts, km$total.iterations, "\n")

set.seed(42)
km <- ewkm(ds, 10, maxrestart=50)
print(km)
cat(km$iterations, km$restarts, km$total.iterations, "\n")

set.seed(42)
km <- ewkm(ds, 10, maxrestart=-1)
print(km)
cat(km$iterations, km$restarts, km$total.iterations, "\n")

set.seed(42)
km <- ewkm(ds, 10, maxiter=1000, maxrestart=0)
print(km)
cat(km$iterations, km$restarts, km$total.iterations, "\n")
