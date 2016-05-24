library(dequer)

n <- 100

l <- list()
for (i in 1:n) l[[i]] <- i

dl <- deque()
for (i in 1:n) pushback(dl, i)

dl.rev <- rev(dl)
l.rev <- rev(l)

stopifnot(all.equal(l.rev, as.list(dl)))

