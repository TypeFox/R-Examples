
setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...) .Call("hiplar_dgeMatrix_norm", x, type, PACKAGE = "HiPLARM"),
	  valueClass = "numeric")

setMethod("norm", signature(x = "matrix", type = "character"),
      function(x, type, ...) .Call("hiplar_dgeMatrix_norm", as(x,"dgeMatrix"), type, PACKAGE = "HiPLARM"),
      valueClass = "numeric")


setMethod("rcond", signature(x = "dgeMatrix", norm = "character"),
	  function(x, norm, ...)  {
	      if({d <- dim(x); d[1] == d[2]})
		  .Call("hiplar_dgeMatrix_rcond", x, norm, PACKAGE = "HiPLARM")
	      else rcond(qr.R(qr(if(d[1] < d[2]) t(x) else x)), norm=norm, ...)
	  },
	  valueClass = "numeric")

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call("hiplar_dgeMatrix_solve", a, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")


setMethod("solve", signature(a = "dgeMatrix", b = "ddenseMatrix"),
	  function(a, b, ...) .Call("hiplar_dgeMatrix_matrix_solve", a, b, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) .Call("hiplar_dgeMatrix_matrix_solve", a, b, PACKAGE = "HiPLARM"),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "sparseMatrix"),
	  function(a, b, ...) .Call("hiplar_dgeMatrix_matrix_solve", a,
				    as(b, "denseMatrix"), PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")


setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, warnSing = TRUE, ...) .Call("hiplar_dgeMatrix_LU", x, warnSing, PACKAGE = "HiPLARM"),
	  valueClass = "denseLU")

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  .Call("hiplar_dgeMatrix_determinant", x, TRUE, PACKAGE = "HiPLARM"))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
	  .Call("hiplar_dgeMatrix_determinant", x, logarithm, PACKAGE = "HiPLARM"))

