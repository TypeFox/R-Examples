# Compute time
make.time <- function(object){
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  paste(h, "h:", m, "m:", s, "s", sep = "")
}

repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

repCols <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

make.add <- function(row, col, Ka){
  sa <- makeL(1:rep(0.5 * Ka * (Ka + 1)))
  for (k in row:col){ 
    cb <- sa[k, row]
    form <- paste(paste("x",  cb:(cb + (Ka - k)), sep = ""), paste("x", cb, sep = ""), sep = "*")
  }
  form
}

#Function repmat like Matlab
repmat <- function(x, dimen){
  kronecker(array(1, dim = dimen), x)
}


