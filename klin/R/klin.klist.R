`klin.klist` <-
function(A) {
  len <- length(A)
  if (len==0)
    Matrix(1,1,1)
  else if (len==1)
    A[[1]]
  else {
    B <- A[[1]]
    for (i in 2:len)
      B <- B %x% A[[i]]
    B
  }
}

