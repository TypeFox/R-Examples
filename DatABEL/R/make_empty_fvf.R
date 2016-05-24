#' makes empty filevector object
#'
#' function to generate empty filevector object (and disk files)
#'
#' @param name name fo the file to be assoiated with new object
#' @param nvariables number of variables (R columns)
#' @param nobservations number of observations (R rows)
#' @param type data type of the object ("UNSIGNED_SHORT_INT",
#' "SHORT_INT", "UNSIGNED_INT", "INT", "FLOAT", "DOUBLE", "CHAR", "UNSIGNED_CHAR")
#' @param cachesizeMb what cache size to use for newly generated 'databel' object
#' @param readonly whether to open new 'databel' in readonly mode
#'
#' @return databel object; also file is created in file system
#' @export
#'

make_empty_fvf <- function(name, nvariables, nobservations,type =
                           "DOUBLE", cachesizeMb = 64, readonly =
                           FALSE)
{
    if (!is(name,"character") || length(name) > 1) {
        stop("file name argument is either non-character or non-scalar");
    }
    nvariables    <- as.integer(nvariables)
    nobservations <- as.integer(nobservations)
    if (nvariables <= 0) stop("number of variables should be > 0 ");
    if (nobservations <= 0) stop("number of observations should be > 0 ");


    intype <- filevector_type(type)

    result <- .Call("ini_empty_FileMatrix_R",
                    as.character(name),
                    as.integer(nvariables),
                    as.integer(nobservations),
                    as.integer(intype));

    if (!result) stop(paste("failed to make filevector file", name))
    # print(as.character(name))
    return(databel(as.character(name),
                   cachesizeMb=cachesizeMb,
                   readonly=readonly))
}
