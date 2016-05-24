left <- function(x0,y0,x1,y1,x2,y2){
  n0 <- length(x0)
  if(length(y0)!=n0)
    stop("x0 and y0 have different length!")
  if(length(x1)!=1 |
     length(y1)!=1 |
     length(x2)!=1 |
     length(x2)!=1)
    stop("x1,y1,x2,y2 should all have length 1!")
     
  
  .Fortran("vleft",
           as.double(x1),
           as.double(y1),
           as.double(x2),
           as.double(y2),
           as.double(x0),
           as.double(y0),
           as.integer(n0),
           left=logical(n0),
           PACKAGE="tripack")$left
}

