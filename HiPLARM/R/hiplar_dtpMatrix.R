
setMethod("solve", signature(a = "dtpMatrix", b="ddenseMatrix"),
	  function(a, b, ...) .Call("magma_dtpMatrix_matrix_solve", a, b, PACKAGE="HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtpMatrix", b="matrix"),
	  function(a, b, ...) .Call("magma_dtpMatrix_matrix_solve", a, b, PACKAGE="HiPLARM"),
	  valueClass = "dgeMatrix")
