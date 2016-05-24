recodeVar <- function(x, src, tgt, default=NULL, keep.na=TRUE){

  if (length(src)!=length(tgt)){
    stop("length of src not equal to length of tgt")
  }
  mtc <- lapply(src, function(zzz){which(x %in% zzz)})
  idx <- seq_along(x)
  unmatch <- setdiff(idx, unlist(mtc))
  
  if (is.factor(x)){
    val <- as.character(x)
  } else {
    val <- x
  }
  for (ii in 1:length(tgt))
    val[mtc[[ii]]] <- tgt[[ii]]

  if (!is.null(default)){
    if (keep.na){
      iii <- intersect(which(!is.na(x)), unmatch)
      val[iii] <- default
    } else {
      val[unmatch] <- default
    }
  }
  
  if (is.factor(x))
    val <- as.factor(val)
  val

  return(val)
}
