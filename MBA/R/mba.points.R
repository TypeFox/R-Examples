"mba.points" <- function(xyz, xy.est, n = 1, m = 1, h = 8, extend=TRUE, verbose=TRUE, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not a parameter")
  }
  
  if(missing(xyz)){stop("error: xyz matrix or data frame must be specified")}
  if(missing(xy.est)){stop("error: xy.est matrix or data frame must be specified")}
  
  if(!any(is.matrix(xyz), is.data.frame(xyz))){stop("error: xyz must be a matrix or data frame")}
  if(any(ncol(xyz) != 3, nrow(xyz) == 0)){stop("error: xyz must have 3 columns corresponding to x, y, z and at least one row")}

  if(!any(is.matrix(xy.est), is.data.frame(xy.est))){stop("error: xy.est must be a matrix or data frame")}
  if(any(ncol(xy.est) != 2, nrow(xy.est) == 0)){stop("error: xy.est must have 2 columns corresponding to x, y coordinates of the point for which to predict z")}

  if(m <= 0){stop("error: m must be a positive integer")}
  if(n <= 0){stop("error: n must be a positive integer")}
  if(h <= 0){stop("error: h must be a positive integer")}

  xyz <- as.matrix(xyz)
  xy.est <- as.matrix(xy.est)
  storage.mode(xyz) <- storage.mode(xy.est) <- "double"
  
  out <- .Call("MBAPoints", xyz, xy.est,  as.integer(m), as.integer(n), as.integer(h), as.integer(extend), as.integer(verbose))
  
  out1 <- list()
  out1$xyz.est <- cbind(xy.est, out)
  colnames(out1$xyz.est) <- c("x","y","z")
  out1
  
}
