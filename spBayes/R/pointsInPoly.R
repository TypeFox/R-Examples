pointsInPoly <- function(poly, points, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(poly)){stop("error: poly must be specified")}
  if(missing(points)){stop("error: points must be specified")}  

  if(!is.matrix(poly)){stop("error: poly must be a matrix")}
  if(ncol(poly) != 2){stop("error: poly must be a matrix with two columns")}

  if(!is.matrix(points)){stop("error: points must be an matrix")}
  if(ncol(points) != 2){stop("error: pointns must be a matrix with two columns")}
    
  n.verts <- nrow(poly)
  n.pts <- nrow(points)
  in.pt.indx <- rep(0, n.pts)
  n.in.pts <- 0
  
  storage.mode(poly) <- "double"
  storage.mode(n.verts) <- "integer"
  storage.mode(points) <- "double"
  storage.mode(n.pts) <- "integer"
  storage.mode(in.pt.indx) <- "integer"
  storage.mode(n.in.pts) <- "integer"

  junk <- .Call("ptsInPoly", poly, n.verts, points, n.pts, in.pt.indx, n.in.pts);
  
  if(n.in.pts == 0){
    NA
  }else{ 
    1+in.pt.indx[1:n.in.pts]##+1 for R indexing
  }
}
