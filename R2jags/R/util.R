fround <- function (x, digits) {
    format (round (x, digits), nsmall=digits)
}
  
pfround <- function (x, digits) {
    print (fround (x, digits), quote=FALSE)
}
 

m <- NULL

repath <- function(x) {
  n.part <- length(x[[1]])
  temp <- rep(NA, n.part)
  for (i in 1:n.part){
    if (nchar(x[[1]][i])> 8){
      temp[i] <- paste(substr(x[[1]][i], 1, 6), "~1/", sep="")
    }
    if (nchar(x[[1]][i])<=8){
      temp[i] <- paste(x[[1]][i],"/", sep="")
    }
  }
  return(temp)
}
    

win2unixdir <- function(windir){
  Dir <- substr(windir, 1, 3)
  tempdir <- substr(windir, 4, 1000)
  tempdir <- gsub(" ", "", tempdir)
  tempdir <- strsplit(tempdir, "/")
  n.part <- length(tempdir[[1]])
  if (n.part>0){
    tempdir <- repath(tempdir)
    path <- ""
    for (i in 1:n.part){
      path <- paste(path, tempdir[i], sep="")
    }
    newpath <- paste(Dir, path, sep="")
  }
  else {
    newpath <- Dir
  }
  return(c(newpath))
}


read.jagsdata <- function(file)
{
  e <- new.env()
  eval(parse(file), e)
  return(as.list(e))
}
