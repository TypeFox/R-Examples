.onLoad <- function(libname, pkgname) {
	.C("hiplar_init", PACKAGE="HiPLARM")

    if ((.Call("hiplarWithMAGMA", PACKAGE="HiPLARM") == TRUE) &&
        (.Call("hiplarWithPLASMA", PACKAGE="HiPLARM") == TRUE)) {
       	checkFile()
    }
	
    if (.Call("hiplarWithMAGMA", PACKAGE="HiPLARM") == FALSE) {
        setMethod("solve", signature(a = "dtpMatrix", b="ddenseMatrix"),
    	  function(a, b, ...) .Call("dtpMatrix_matrix_solve", a, b),
    	  valueClass = "dgeMatrix")

        setMethod("solve", signature(a = "dtpMatrix", b="matrix"),
    	  function(a, b, ...) .Call("dtpMatrix_matrix_solve", a, b),
    	  valueClass = "dgeMatrix")

        setMethod("%*%", signature(x = "dspMatrix", y = "ddenseMatrix"),
          function(x, y) .Call("dspMatrix_matrix_mm", x, y),
          valueClass = "dgeMatrix")

        setMethod("%*%", signature(x = "dspMatrix", y = "matrix"),
          function(x, y) .Call("dspMatrix_matrix_mm", x, y),
          valueClass = "dgeMatrix")

        setMethod("%*%", signature(x = "dtpMatrix", y = "ddenseMatrix"),
    	  function(x, y) .Call("dtpMatrix_matrix_mm", x, y))

        setMethod("%*%", signature(x = "dgeMatrix", y = "dtpMatrix"),
    	  function(x, y) .Call("dgeMatrix_dtpMatrix_mm", x, y))
    }
}

.onUnload <- function(libname) {
	.C("hiplar_deinit", PACKAGE="HiPLARM")
}

