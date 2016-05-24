`pimg.list` <-
function(times, xlim, ylim, img.dim, Z = TRUE) {

  #epoch <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "GMT")
  #tm <- 12*3600*(round(unclass(times)/(12*3600))) + epoch
  tm <- times

  ## if these are Zs, then the times specify the X times
  n <- length(tm) - Z

  lst <- vector("list",length=n)
  for(i in 1:n) lst[[i]] <- pimg(xlim[1], xlim[2], img.dim[1], ylim[1], ylim[2], img.dim[2])
  
  if (Z) {
  	Ztimes <- tm[-length(tm)] + diff(unclass(tm))/2
  	attr(lst, "times") <- Ztimes
  	attr(lst, "Xtimes") <- tm
  } else {
  	attr(lst, "times") <- tm
  }
  attr(lst, "Z") <- Z
   
  class(lst) <- c("pimg.list", "list")
  lst
}

