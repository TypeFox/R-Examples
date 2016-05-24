### Checking for memory leaks
library(dequer)

n <- 1e5

l1 <- lapply(n:1, identity)

system.time({
  dl <- deque()
  for (i in 1:n) push(dl, i)
  l2 <- as.list(dl)
})

stopifnot(all.equal(l1, l2))

rm(dl, l2)
invisible(gc())

