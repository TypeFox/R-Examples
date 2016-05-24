#' generates temporary file name
#'
#' function to generate temporary file name
#'
#' @param path path to directory where the temporary file will be
#' located
#'
#' @param withFVext whether function should check presence of *FVD and
#' *FVI files too
#' @export
#'

get_temporary_file_name <- function(path=".", withFVext=TRUE)
{
    tmpname <- paste(path, "/tmp",
                     round(runif(1, min=10000, max=1000000)),
                     sep="")

    while (file.exists(tmpname)
           || file.exists(paste(tmpname, ".fvi", sep=""))
           || file.exists(paste(tmpname, ".fvd", sep=""))) {
        tmpname <- paste(path, "/tmp",
                         round(runif(1, min=10000, max=1000000)),
                         sep="")
    }
    return(tmpname)
}
