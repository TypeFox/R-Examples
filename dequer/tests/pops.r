### Checking for memory leaks
library(dequer)

dl <- deque()

push(dl, 1234)
pop(dl)

push(dl, 1)
push(dl, 2)
pushback(dl, 3)

pop(dl)
popback(dl)

l2 <- as.list(dl)

l1 <- list(1)
stopifnot(all.equal(l1, l2))

invisible({rm(l2);gc()})
