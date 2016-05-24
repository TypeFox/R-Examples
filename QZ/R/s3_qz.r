### S3 functions for generalized eigen values.

qz <- function(A, B = NULL, select = NULL, only.values = FALSE, ...){
  if(!is.null(select)){
    only.values <- FALSE
  }

  if(!is.null(B)){
    if(is.complex(A) && is.complex(B)){
      if(only.values){
        ret <- qz.zgges(A, B, vsl = FALSE, vsr = FALSE, ...)
      } else{
        ret <- qz.zgges(A, B, vsl = TRUE, vsr = TRUE, ...)
      }

      if(!is.null(select)){
        ret <- qz.ztgsen(ret$S, ret$T, ret$Q, ret$Z, select, ijob = 4L,
                         want.Q = TRUE, want.Z = TRUE, ...)
      }
    } else if(is.double(A) && is.double(B)){
      if(only.values){
        ret <- qz.dgges(A, B, vsl = FALSE, vsr = FALSE, ...)
      } else{
        ret <- qz.dgges(A, B, vsl = TRUE, vsr = TRUE, ...)
      }

      if(!is.null(select)){
        ret <- qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z, select, ijob = 4L,
                         want.Q = TRUE, want.Z = TRUE, ...)
      }
    } else{
      stop("A and B are not of the same type (either complex or double).")
    }
  } else{
    if(is.complex(A)){
      if(only.values){
        ret <- qz.zgees(A, vs = FALSE, ...)
      } else{
        ret <- qz.zgees(A, vs = TRUE, ...)
      }

      if(!is.null(select)){
        ret <- qz.ztrsen(ret$T, ret$Q, select, job = "B", want.Q = TRUE, ...)
      }
    } else if(is.double(A)){
      if(only.values){
        ret <- qz.dgees(A, vs = FALSE, ...)
      } else{
        ret <- qz.dgees(A, vs = TRUE, ...)
      }

      if(!is.null(select)){
        ret <- qz.dtrsen(ret$T, ret$Q, select, job = "B", want.Q = TRUE, ...)
      }
    } else{
      stop("A is neither complex nor double.")
    }
  }

  ret
} # End of qz().

