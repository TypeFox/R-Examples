print.core <- function(x, vol=NULL, ...) {

  cv <- x$volume; cv[is.na(cv)]=min(cv,na.rm=TRUE)
  cl <- x$length; cl[is.na(cl)]=min(cl,na.rm=TRUE)
  ca <- NULL; cx <- NULL; cr <- NULL

  if(is.null(vol))
    vol <- 1

  ind <- (cv<=vol)
  ca <- sort(x$step.inds[ind])
  cx <- atom2xyz(ca)
  cr <- sort(x$resno[ind])
  nc <- length(ca)
  
  cat(paste("#",nc,
            "positions (cumulative volume <=",
            vol,"Angstrom^3)"),"\n")
  if(nc==0) {
    cat(paste("# Min volume is",round(min(cv),3)),"\n")
  } else {
    print(bounds(as.numeric(cr)), ...)  
  }
  ##NextMethod("print", x, quote = FALSE, right = TRUE, ...)
  invisible(list(atom=ca, xyz=cx, resno=cr))
}
