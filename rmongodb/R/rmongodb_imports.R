#' @useDynLib rmongodb
#' @import jsonlite
#' @import plyr
#'
.onUnload <- function(libpath)
    library.dynam.unload("rmongodb", libpath)


