library(dequer)

f <- function(n, k)
{
  dl1 <- deque()
  for (i in 1:n) pushback(dl1, i)
  
  dl2 <- sep(dl1, k)
  
  l1 <- as.list(dl1)
  l2 <- as.list(dl2)
  
  truth1 <- lapply(1:k, identity)
  truth2 <- lapply((k+1):n, identity)
  
  stopifnot(all.equal(l1, truth1))
  stopifnot(all.equal(l2, truth2))
  
  invisible({rm(dl1, dl2);gc()})
}


f(7, 3)
f(7, 4)
