plotsurvivors <-
function(x,status=rep(1,length(x))) {
  if (any(x==0)) {
    x <- x[-which(x==0)]
    status <- status[-which(x==0)]
  }
  if (any(status==0)) {
    x <- x[-which(status==0)]
    status <- status[-which(status==0)]
  }
  n <- length(x)
  tri <- sort(x)
  tri2 <- unique(tri)
  alive <- integer(length(tri2))
  alive[1] <- n
  for (i in 2:length(unique(tri2))) {
    alive[i] <- alive[i-1]-length(which(tri==tri2[i-1]))
  }
  plot(tri2,alive,pch=16,xlab="Time",ylab="Survivors (log scale)",log="y")
  result <- list(n=n,time=tri2,alive=alive)
  invisible(result)
}
