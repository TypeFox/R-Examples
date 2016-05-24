library(dequer)


n <- 2e5

#cat("-------- push --------\n")

#system.time({
  #l <- list()
  #for (i in 1:n) l[[i]] <- n-i+1
#})


#system.time({
  #dl <- deque()
  #for (i in 1:n) push(dl, i)
  #l2 <- as.list(dl)
#})

#all.equal(l, l2)


#rm(l, dl, l2)
#invisible(gc())



cat("-------- pushback --------\n")

system.time({
  l <- list()
  for (i in 1:n) l[[i]] <- i
})


system.time({
  dl <- deque()
  for (i in 1:n) pushback(dl, i)
  l2 <- as.list(dl)
})

all.equal(l, l2)

rev(dl)
