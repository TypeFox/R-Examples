print.marked.version <- function()
{ library(help=marked)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version))
	{
		um <- strsplit(version," ")[[1]]
	    version <- um[nchar(um)>0][2]
	    hello <- paste("This is marked ",version,"\n",sep="")
	} else
		hello <- "This is marked"
	packageStartupMessage(hello)
}
.onAttach <- function(...) { 
	print.marked.version()
}

.onUnload <- function(libpath)
{
#  unloadNamespace("marked")
#  library.dynam.unload("marked", libpath)
  cat("\nBye-Bye from marked\n\n")
  return(invisible())
}
#' @import methods
### setMethods for Matrix package
setMethod(f=crossprod, signature=signature(x="dgeMatrix", y="dtTMatrix"), function(x,y) callNextMethod())

# setMethod(f="%*%", signature=signature(x="ddiMatrix", y="dtTMatrix"), function(x,y) Matrix:::diagCspprod(as(x, "CsparseMatrix"), y))
# setMethod(f=crossprod, signature=signature(x="dgeMatrix", y="dtTMatrix"), 
#           definition=function (x, y = NULL){
#             t(.Call("Csparse_dense_crossprod", as(y, "CsparseMatrix"), x, PACKAGE="Matrix"))
#           }
# )
# 

