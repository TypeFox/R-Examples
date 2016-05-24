library(dequer)

dl <- deque()
push(dl, 1234)
push(dl, 897)
pushback(dl, "asdf")

stopifnot(length(dl) == 3)

l1 <- list(897, 1234, "asdf")
l2 <- as.list(dl)

stopifnot(all.equal(l1, l2))

