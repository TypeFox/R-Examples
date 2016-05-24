co_table <- function(x, margin, collapse = ".")
{
  if (!is.array(x)) 
      stop("x is not an array")
  if("xtabs" %in% class(x)) attr(x, "call") <- NULL

  dx <- dim(x)
  idx <- lapply(dx, function(i) 1:i)
  dn <- dimnames(x)
  
  if(is.character(margin)) {
    if(is.null(dn)) stop("margin must be an index when no dimnames are given")
    margin <- which(names(dn) %in% margin)
  }

  idxm <- expand.grid(idx[margin])    
  cotab1 <- function(i) {
    idx[margin] <- lapply(1:length(margin), function(j) idxm[i,j])
    rval <- as.table(do.call("[", c(list(x), idx, list(drop = FALSE))))
    if(length(dim(rval)) > 1) {
      dim(rval) <- dim(x)[-margin]
      dimnames(rval) <- dimnames(x)[-margin]
    }
    return(rval)
  }    
  rval <- lapply(1:NROW(idxm), cotab1)
  if(!is.null(dn)) names(rval) <- apply(expand.grid(dn[margin]), 1, function(z) paste(z, collapse = collapse))    

  return(rval)
}
