
setMethod("chol", signature(x = "dpoMatrix"),
	  function(x, pivot, ...) .Call("hiplar_dpoMatrix_chol", x, PACKAGE = "HiPLARM"))

setMethod("rcond", signature(x = "dpoMatrix", norm = "character"),
          function(x, norm, ...) .Call("hiplar_dpoMatrix_rcond", x, norm, PACKAGE = "HiPLARM"))

setMethod("rcond", signature(x = "dpoMatrix", norm = "missing"),
          function(x, norm, ...) .Call("hiplar_dpoMatrix_rcond", x, "O", PACKAGE = "HiPLARM"))

setMethod("solve", signature(a = "dpoMatrix", b = "missing"),
          function(a, b, ...) .Call("hiplar_dpoMatrix_solve", a, PACKAGE = "HiPLARM"),
          valueClass = "dpoMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "dgeMatrix"),
          function(a, b, ...) .Call("hiplar_dpoMatrix_dgeMatrix_solve", a, b, PACKAGE = "HiPLARM"),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "matrix"),
          function(a, b, ...) .Call("hiplar_dpoMatrix_matrix_solve", a, b, PACKAGE = "HiPLARM"),
          valueClass = "matrix")

