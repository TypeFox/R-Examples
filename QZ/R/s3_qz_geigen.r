### Generalized eigen values.

qz.geigen <- function(A, B = NULL, only.values = FALSE, ...){
  if(!is.null(B)){
    if(is.complex(A) && is.complex(B)){
      if(only.values){
        ret <- qz.zggev(A, B, vl = FALSE, vr = FALSE, ...)
      } else{
        ret <- qz.zggev(A, B, vl = TRUE, vr = TRUE, ...)
      }
    } else if(is.double(A) && is.double(B)){
      if(only.values){
        ret <- qz.dggev(A, B, vl = FALSE, vr = FALSE, ...)
      } else{
        ret <- qz.dggev(A, B, vl = TRUE, vr = TRUE, ...)
      }
    } else{
      stop("A and B are not of the same type (either complex or double).")
    }
  } else{
    if(is.complex(A)){
      if(only.values){
        ret <- qz.zgeev(A, vl = FALSE, vr = FALSE, ...)
      } else{
        ret <- qz.zgeev(A, vl = TRUE, vr = TRUE, ...)
      }
    } else if(is.double(A)){
      if(only.values){
        ret <- qz.dgeev(A, vl = FALSE, vr = FALSE, ...)
      } else{
        ret <- qz.dgeev(A, vl = TRUE, vr = TRUE, ...)
      }
    } else{
      stop("A is neither complex nor double.")
    }
  }

  ret
} # End of qz.geigen().

geigen <- qz.geigen
