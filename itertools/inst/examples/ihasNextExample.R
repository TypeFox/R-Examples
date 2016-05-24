library(itertools)

n <- 100
it <- ihasNext(icount(n))

total <- 0
while (hasNext(it))
  total <- total + nextElem(it)

print(total == sum(seq(length=n)))
