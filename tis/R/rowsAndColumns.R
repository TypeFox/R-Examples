columns <- function(z){
  ## returns a list of the columns in z
  x <- as.matrix(z)
  n <- ncol(x)
  xlist <- vector("list", n)
  for(i in 1:n)
    xlist[[i]] <- x[,i]
  if(length(cn <- colnames(x)) == n)
    names(xlist) <- cn
  xlist
}

rows <- function(z){
  ## returns a list of the rows in z
  x <- as.matrix(z)
  n <- nrow(x)
  xlist <- vector("list", n)
  for(i in 1:n)
    xlist[[i]] <- x[i,]
  if(length(rn <- rownames(x)) == n)
    names(xlist) <- rn
  xlist
}
