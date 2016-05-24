library(pcalg)

my.stop <- FALSE
tmp.set <- c(1,2,3)
res <- NULL

while(!my.stop) {
  tmp <- getNextSet(5,3,tmp.set)
  tmp.set <- tmp$nextSet
  if (tmp$wasLast) {
    my.stop <- TRUE
  } else {
    res <- rbind(res,tmp.set)
  }
}

resTrue <- rbind(c(1,2,4),c(1,2,5),c(1,3,4),c(1,3,5),c(1,4,5),c(2,3,4),c(2,3,5),c(2,4,5),c(3,4,5))

if(any(resTrue!=res)) {
  stop("Test of getNextSet: Theoretical values not reproduced!")
}
