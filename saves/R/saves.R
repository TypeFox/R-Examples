#' Save the variables of a data.frame in distinct binary files
#'
#' \code{saves} does what the name suggests: it saves dataframe(s) or list(s) to disk in a special, binary format. This binary format consists of distinct binary files of all separate variables of a dataframe/list merged into an uncompressed tar archive. This is done via a loop, which saves each variable/column to an an external representation of the R objects via \code{save} in a temporary directory. Theese 'RData' files are archived to an 'RDatas' tar archive, uncompressed for better speed.
#' @param ... R objects: the names of the objects to be saved (as symbols or character strings)
#' @param list character vector: the name(s) of the data frame(s) or list(s) to save
#' @param file character vector: the (RDatas) filename(s) in which to save the variables in the current working directory
#' @param overwrite boolean: if TRUE, existing files will be deleted before saving. Default set to FALSE, which will report error on conflicting file names.
#' @param ultra.fast boolean: if TRUE, ultra fast (...) processing is done without any check to parameters, also no archiving or compression is done. Be sure if using this setting, as many uncompressed files could be generated in the working directory's subdirectory named to \code{df}. Only recommended for servers dealing with lot of R objects' saves and loads in a monitored environment.
#' @return The saved filename(s) (invisible).
#' @export
#' @seealso \code{loads} to load R objects from RDatas binary format
#' @examples \dontrun{
#' ## Saving the demo dataset to evs.2000.hun.RDatas in current working directory.
#' data(evs.2000.hun)
#' saves(evs.2000.hun)
#' ## Saving both the demo dataset and mtcars to current working directory
#' saves(evs.2000.hun, mtcars)
#' saves(list=c('evs.2000.hun', 'mtcars'))
#' ## Saving all kind of cars :)
#' saves(cars, mtcars, overwrite = T)
#' saves(list=c('cars', 'mtcars'), overwrite = T)
#' }
saves <- function (..., list=character(), file=NULL, overwrite=FALSE, ultra.fast=FALSE) {

    names <- as.character(substitute(list(...)))[-1L]
    list <- c(list, names)

    if (ultra.fast == TRUE) {
        df <- list[1]
        data <- get(df)
        dir.create(df)
        e <- as.environment(data)
        for (i in 1:length(data)) {
            save(list=names(data)[i], file=paste(df, '/', names(data)[i], '.RData', sep=''), compress=FALSE, precheck=FALSE, envir = e)
        }
        return(invisible(df))
    }

    if (is.null(file)) file <- paste(list, '.RDatas', sep='')
    if (length(list) != length(file)) stop('Bad number of files given!')

    for (i in 1:length(list)) {
        if (inherits(try(data <- get(list[i]), silent=TRUE), "try-error")) stop(paste('No dataframe/list given or `', list[i] ,'` is not a dataframe/list!'))
        if (file.exists(file[i])) {
            if (overwrite == TRUE) {
                file.remove(file[i])
            } else {
                stop(paste('Destination filename `', file[i], '` already exists! Use other filename or use paramater `overwrite` set to TRUE.'))
            }
        }
        if (!is.data.frame(data) & !is.list(data)) stop(paste('No dataframe/list given or `', list[i] ,'` is not a dataframe/list!'))

        tmp <- tempfile('saves.dir-')
        dir.create(tmp)
        e <- as.environment(data)
        lapply(names(data), function(x) save(list=x, file=paste(tmp, '/', x, '.RData', sep=''), envir = e))
        w <- getwd()
        setwd(tmp)
        tar(paste(w, '/', file[i], sep=''), '.', compression='none')
        setwd(w)
        unlink(tmp, recursive = TRUE)
    }
    invisible(file)
}
