setMethod("norm", signature(x = "ldenseMatrix", type = "character"),
	  function(x, type, ...)
          .Call("hiplar_dgeMatrix_norm", as(as(x,"dMatrix"),"dgeMatrix"), type, PACKAGE = "HiPLARM"),
	  valueClass = "numeric")
