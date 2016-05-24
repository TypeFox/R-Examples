#' Generic function 
#' @name createSimulationTable
#' @aliases createSimulationTable
#' @title Generic function 
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{createSimulationTable.ezsim}} 
createSimulationTable <-
function(x,...){
    UseMethod("createSimulationTable",x)
}

#' Generic function
#' @name generate
#' @aliases generate
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{generate.parameterDef}}
generate <-
function(x,...){
    UseMethod("generate",x)
}

#' Generic function
#' @name setBanker
#' @aliases setBanker
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{setBanker.parameterDef}}
setBanker <-
function(x,...){
    UseMethod("setBanker",x)
}

#' Generic function
#' @name run
#' @aliases run
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{run.ezsim}}

run <-
function(x,...){
    UseMethod("run",x)
}
#' Generic function
#' @name getSelectionName
#' @aliases getSelectionName
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{getSelectionName.ezsim}}, \code{\link{getSelectionName.summary.ezsim}}

getSelectionName <-
function(x,...){
    UseMethod("getSelectionName",x)
}
#' Generic function
#' @name setBanker
#' @aliases setBanker
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @seealso \code{\link{setBanker.parameterDef}}

setBanker <-
function(x,...){
    UseMethod("setBanker",x)
}

#' Generic function
#' @name setSelection
#' @aliases setSelection
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @keywords internal
#' @export
#' @seealso \code{\link{setSelection.parameterDef}}
setSelection <-
function(x,...){
    UseMethod("setSelection",x)
}

#' Generic function
#' @name test
#' @aliases test
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @keywords internal
#' @export
#' @seealso \code{\link{test.ezsim}}
test <-
function(x,...){
    UseMethod("test",x)
}

#' Generic function
#' @name createSimulationTable
#' @aliases createSimulationTable
#' @title Generic function
#' @param x Object
#' @param \dots Further arguments
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @keywords internal
#' @export
createSimulationTable <-
function(x,...){
    UseMethod("createSimulationTable",x)
}

