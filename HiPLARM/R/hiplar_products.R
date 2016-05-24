
## General method for dense matrix multiplication in case specific methods
## have not been defined.
setMethod("%*%", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call("hiplar_dgeMatrix_matrix_mm",
			       .Call("dup_mMatrix_as_dgeMatrix", x), y, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("hiplar_dgeMatrix_matrix_mm", x, y, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) .Call("hiplar_dgeMatrix_matrix_mm", x, y, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) .Call("hiplar_dgeMatrix_matrix_mm", y, x, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

.dsy_m_mm <- function(x, y) .Call("hiplar_dsyMatrix_matrix_mm", x, y, FALSE, PACKAGE = "HiPLARM")
setMethod("%*%", signature(x = "dsyMatrix", y = "matrix"),  .dsy_m_mm)
setMethod("%*%", signature(x = "dsyMatrix", y = "ddenseMatrix"),  .dsy_m_mm)
## for disambiguity :
setMethod("%*%", signature(x = "dsyMatrix", y = "dsyMatrix"),  .dsy_m_mm)

setMethod("%*%", signature(x = "ddenseMatrix", y = "dsyMatrix"),
          function(x, y) .Call("hiplar_dsyMatrix_matrix_mm", y, x, TRUE, PACKAGE = "HiPLARM"))
setMethod("%*%", signature(x = "matrix", y = "dsyMatrix"),
          function(x, y) .Call("hiplar_dsyMatrix_matrix_mm", y, x, TRUE, PACKAGE = "HiPLARM"))


setMethod("%*%", signature(x = "dspMatrix", y = "ddenseMatrix"),
          function(x, y) .Call("magma_dspMatrix_matrix_mm", x, y, PACKAGE = "HiPLARM"),
          valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "dspMatrix", y = "matrix"),
          function(x, y) .Call("magma_dspMatrix_matrix_mm", x, y, PACKAGE = "HiPLARM"),
          valueClass = "dgeMatrix")


setMethod("%*%", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_dtrMatrix_mm", x, y, FALSE, FALSE, PACKAGE = "HiPLARM"))

setMethod("%*%", signature(x = "dtrMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", x, y, FALSE, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", x, y, FALSE, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "ddenseMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", y, x, TRUE, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dtrMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", y, x, TRUE, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")


setMethod("%*%", signature(x = "dtpMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call("magma_dtpMatrix_matrix_mm", x, y, PACKAGE = "HiPLARM"))
setMethod("%*%", signature(x = "dgeMatrix", y = "dtpMatrix"),
	  function(x, y) .Call("magma_dgeMatrix_dtpMatrix_mm", x, y, PACKAGE = "HiPLARM"))



###--- II --- crossprod -----------------------------------------------------

setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_crossprod", x, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dpoMatrix")

## crossprod (x,y)
setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_dgeMatrix_crossprod", x, y, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_matrix_crossprod", x, y, FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call("hiplar_dgeMatrix_matrix_crossprod", x, as.matrix(as.double(y)), FALSE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

## "dtrMatrix" - remaining (uni)triangular if possible
setMethod("crossprod", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_dtrMatrix_mm", x, y, FALSE, TRUE, PACKAGE = "HiPLARM"))

setMethod("crossprod", signature(x = "dtrMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", x, y, FALSE, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", x, y, FALSE, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")


###--- III --- tcrossprod ---------------------------------------------------

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_dgeMatrix_crossprod", x, y, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_matrix_crossprod", x, y, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call("hiplar_dgeMatrix_matrix_crossprod", x, rbind(as.double(y)), TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dgeMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call("hiplar_dgeMatrix_crossprod", x, TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dpoMatrix")

if(FALSE) { ## this would mask 'base::tcrossprod'
setMethod("tcrossprod", signature(x = "matrix", y = "missing"),
	  function(x, y = NULL)
	  .Call("hiplar_dgeMatrix_crossprod", as(x, "dgeMatrix"), TRUE, PACKAGE = "HiPLARM"),
	  valueClass = "dpoMatrix")
}


setMethod("tcrossprod", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("hiplar_dtrMatrix_dtrMatrix_mm", y, x, TRUE, TRUE, PACKAGE = "HiPLARM"))

setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "dtrMatrix"),
 	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", y, x, TRUE, TRUE, PACKAGE = "HiPLARM"))

setMethod("tcrossprod", signature(x = "matrix", y = "dtrMatrix"),
 	  function(x, y) .Call("hiplar_dtrMatrix_matrix_mm", y, x, TRUE, TRUE, PACKAGE = "HiPLARM"))

