library(dequer)

n <- 3
m <- 4

dl <- deque()
for (i in 1:n) push(dl, i)

dl2 <- deque()
for (i in 1:m) push(dl2, i)

combine(dl, dl2)

l <- as.list(dl)

truth <- list(3, 2, 1, 4, 3, 2, 1)

stopifnot(all.equal(l, truth))

invisible({rm(dl, dl2);gc()})
