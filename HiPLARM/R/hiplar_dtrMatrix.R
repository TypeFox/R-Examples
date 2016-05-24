#
#setMethod("determinant", signature(x = "dtrMatrix", logarithm = "missing"),
#	  function(x, logarithm, ...) callGeneric(x, TRUE))
#
#setMethod("determinant", signature(x = "dtrMatrix", logarithm = "logical"),
#	  function(x, logarithm, ...) mkDet(diag(x), logarithm))
#
#setMethod("norm", signature(x = "dtrMatrix", type = "character"),
#	  function(x, type, ...)
#	  .Call(dtrMatrix_norm, x, type),
#	  valueClass = "numeric")
#
#setMethod("norm", signature(x = "dtrMatrix", type = "missing"),
#	  function(x, type, ...)
#	  .Call(dtrMatrix_norm, x, "O"),
#	  valueClass = "numeric")
#
#setMethod("rcond", signature(x = "dtrMatrix", norm = "character"),
#	  function(x, norm, ...)
#	  .Call(dtrMatrix_rcond, x, norm),
#	  valueClass = "numeric")
#
#setMethod("rcond", signature(x = "dtrMatrix", norm = "missing"),
#	  function(x, norm, ...)
#	  .Call(dtrMatrix_rcond, x, "O"),
#	  valueClass = "numeric")

setMethod("chol2inv", signature(x = "dtrMatrix"),
	  function (x, ...) {
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      if (x@diag != "N") x <- diagU2N(x)
	      .Call("hiplar_dtrMatrix_chol2inv", x, PACKAGE = "HiPLARM")
	  })

setMethod("solve", signature(a = "dtrMatrix", b="missing"),
	  function(a, b, ...) {
	      ## warn, as e.g. CHMfactor have 'system' as third argument
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      .Call("hiplar_dtrMatrix_solve", a, PACKAGE = "HiPLARM")
	  }, valueClass = "dtrMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="ddenseMatrix"),
	  function(a, b, ...) {
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      .Call("hiplar_dtrMatrix_matrix_solve", a, b, PACKAGE = "HiPLARM")
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="dMatrix"),
	  function(a, b, ...) {
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      .Call("hiplar_dtrMatrix_matrix_solve", a, as(b,"denseMatrix"), PACKAGE = "HiPLARM")
	  }, valueClass = "dgeMatrix")
setMethod("solve", signature(a = "dtrMatrix", b="Matrix"),
	  function(a, b, ...) {
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      .Call("hiplar_dtrMatrix_matrix_solve", a, as(as(b, "dMatrix"),
						  "denseMatrix"), PACKAGE = "HiPLARM")
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="matrix"),
	  function(a, b, ...) {
	      if(length(list(...)))
		  warning("arguments in ",deparse(list(...))," are disregarded")
	      .Call("hiplar_dtrMatrix_matrix_solve", a, b, PACKAGE = "HiPLARM")
	  }, valueClass = "dgeMatrix")

