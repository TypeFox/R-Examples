outer.convhull<-function(cx,cy,px,py,FUN,duplicate="remove",...)
  {
    nx<-length(cx)
    ny<-length(cy)
    np<-length(px)
    if(length(py)!=np)
      stop("length of cx and cy differ")
    if (is.character(FUN)) 
      FUN <- get(FUN, mode = "function", inherits = TRUE)
    p.tr<-tri.mesh(px,py,duplicate)

    ans<-matrix(FUN(matrix(cx, nx, ny),
                    matrix(cy, nx, ny, byrow = TRUE), 
                    ...), nx, ny)
    ans[!in.convex.hull(p.tr,matrix(cx, nx, ny),
                        matrix(cy, nx, ny, byrow = TRUE))]<-NA
    ans
  }
