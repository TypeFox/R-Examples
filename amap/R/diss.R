diss <- function (x, w=rep(1,ncol(x)))
  {
    n <- nrow(x)
    p <- ncol(x)
    if(length(w) != p)
      {
        warning("Error in dimention on either w or x")
        return(NULL)
      }
    
    res <- .C("diss",
              as.integer(x),
              double(n*n),
              n,p,
              as.double(w),
              PACKAGE="amap")

    matrix(res[[2]],n)

  }
