#
#setMethod("rcond", signature(x = "dsyMatrix", norm = "character"),
#          function(x, norm, ...)
#          .Call(dsyMatrix_rcond, x, norm),
#          valueClass = "numeric")
#
#setMethod("rcond", signature(x = "dsyMatrix", norm = "missing"),
#          function(x, norm, ...)
#          .Call(dsyMatrix_rcond, x, "O"),
#          valueClass = "numeric")
#
#setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
#          function(a, b, ...) .Call(dsyMatrix_solve, a),
#          valueClass = "dsyMatrix")
#
#setMethod("solve", signature(a = "dsyMatrix", b = "matrix"),
#          function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b),
#          valueClass = "dgeMatrix")
#
#setMethod("solve", signature(a = "dsyMatrix", b = "ddenseMatrix"),
#	  function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b))
#setMethod("solve", signature(a = "dsyMatrix", b = "denseMatrix"), ## eg. for ddi* or ldi*
#	  function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, as(b,"dMatrix")))

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...) .Call("hiplar_dsyMatrix_norm", x, type, PACKAGE = "HiPLARM"),
          valueClass = "numeric")

setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call("hiplar_dsyMatrix_norm", x, "O", PACKAGE = "HiPLARM"),
          valueClass = "numeric")

#setMethod("BunchKaufman", signature(x = "dsyMatrix"),
#	  function(x) .Call(dsyMatrix_trf, x))
#

setAs("dsyMatrix", "dpoMatrix",
      function(from){
      if(is.null(tryCatch(.Call("hiplar_dpoMatrix_chol", from, PACKAGE = "HiPLARM"),
                  error = function(e) NULL))) 
          stop("not a positive definite matrix")
      ## else 
      Matrix:::copyClass(from, "dpoMatrix",
            sNames = c("x", "Dim", "Dimnames", "uplo", "factors"))
      })      

