#' Read Fixed Width Format Files containing lines of different Type
#' 
#' Read a table of fixed width formatted data of different types into a \link{data.frame} for each type.
#' 
#' @param file the name of the file which the data are to be read from.
#' 
#' Alternatively, file can be a \link{connection}, which will be opened if necessary, and if so closed at the end of the function call.
#' @param multi.specs A named list of data.frames containing the following columns:
#'     \tabular{ll}{
#'         widths \tab   see \link{read.fwf}\cr
#'         col.names \tab see \link{read.table}\cr
#'         row.names \tab see \link{read.table}\cr
#'     }
#' 
#' For more info on these fields see \code{\link{read.fwf}}.
#' 
#' Note that each list item should have a name. This is important for the select function.
#' @param select A function to select the type of a line. This selector should have parameters:
#'     \tabular{ll}{
#'         line \tab   the line\cr
#'         specs \tab  the multi.specs list that was passed to read.multi.fwf\cr
#'     }
#'     
#' The select function should return the name of the spec that matches the line. read.multi.fwf will then use this name to select the a spec from the passed multi.spec. This is why multi.spec should be a named list. 
#' If there is no match then NULL can be returned.
#' @param header a logical value indicating whether the file contains the names of the variables as its first line. If present, the names must be delimited by sep.
#' @param sep character; the separator used internally; should be a character that does not occur in the file (except in the header).
#' @param skip number of initial lines to skip; see \link{read.fwf}.
#' @param n the maximum number of records (lines) to be read, defaulting to no limit.
#' @param buffersize Maximum number of lines to read at one time
#' @param ... further arguments to be passed to \link{read.fwf}.
#' @return Return value is a named list with an item for each spec in multi.spec. If there was at least one line in file, matching a spec, then the named item will be a \link{data.frame}. Otherwise it will be NULL.
#' @details
#' Known bugs:
#' Warnings on connections that are left open. Haven't figured this out yet. Somehow some files are left opened.
#' @author
#' Panos Rontogiannis \email{p.g.ronto@@gmail.com}
#' @seealso \code{\link{read.fwf}}
#' @examples
#' ff <- tempfile()
#' cat(file = ff, '123456', '287654', '198765', sep = '\n')
#' specs <- list()
#' specs[['sp1']] = data.frame(widths = c(1, 2, 3), 
#'                             col.names = c('Col1', 'Col2', 'Col3'))
#' specs[['sp2']] = data.frame(widths = c(3, 2, 1), 
#'                             col.names = c('C1', 'C2', 'C3'))
#' 
#' myselector <- function(line, specs) {
#'     s <- substr(line, 1, 1)
#'     spec_name = ''
#'     if (s == '1')
#'         spec_name = 'sp1'
#'     else if (s == '2')
#'         spec_name = 'sp2'
#' 
#'     spec_name
#' }
#' 
#' read.multi.fwf(ff, multi.specs = specs, select = myselector)    
#' #> sp1: 1 23 456 \ 1 98 765, sp2: 287 65 4
#' 
#' unlink(ff)
#' @export
read.multi.fwf <- function(file,
                           multi.specs,
                           select,
                           header = FALSE,
                           sep = "\t",
                           skip = 0,
                           n = -1,
                           buffersize = 2000,
                           ...) {
    ##### Helper functions #############################

    # Convert each spec into a list containing a SPEC field. Also add
    # the FILENAME field for the temp file where all lines of this type (spec)
    # will be stored (temporarily) for read.fwf.
    prepareAsList <- function(s) {
        s <- list(SPEC = s)
        s$FILENAME <- tempfile("Rmfwf.")
        return(s)
    }

    # Add header line to a temp file.
    addHeader <- function(s, headerline) {
        cat(file = s$FILENAME, headerline, "\n")
        return(s)
    }

    # Parse a line and write to the temp file of the matching spec.
    doone <- function(line) {
        spec_name <- select(line, multi.specs)
        s <- extended.specs[[spec_name]]
        if (is.null(s))
            return()
        cat(file = s$FILENAME, line, "\n", append = TRUE)
        invisible()
    }

    # Read Fixed-Width Format temp file.
    readFwf <- function(s) {
        fi <- file.info(s$FILENAME)
        if (!is.na(fi$size) & fi$size > 0) {
            s$Data <- read.fwf(file = s$FILENAME,
                               widths = s$SPEC$widths,
                               header = header,
                               sep = sep,
                               row.names = s$SPEC$row.names,
                               col.names = s$SPEC$col.names,
                               ...)
        }
        else {
            #s$Data <- NA
            s <- NULL
        }
        return(s)
    }

    # Prepare an element for return.
    prepareRetval <- function(s) {
        s <- s$Data
        return(s)
    }
    ####################################################

    extended.specs <- lapply(multi.specs, FUN = prepareAsList)

    if (is.character(file)) {
        file <- file(file, "rt")
        on.exit(close(file), add = TRUE)
    }
    else if (!isOpen(file)) {
        open(file, "rt")
        on.exit(close(file), add = TRUE)
    }

    # Handle header line.
    if (header) {
        headerline <- readLines(file, n = 1L)
        lapply(extended.specs, FUN = addHeader, headerline)
    }
    
    # Skip lines.
    if (skip)
        readLines(file, n = skip)

    repeat ({
        if (n == 0L)
            break
        if (n == -1L)
            n <- 16

        raw <- readLines(file, n = n)
        nread <- length(raw)
        if (nread == 0)
            break

        lapply(raw, FUN = doone)
    })

    extended.specs <- lapply(extended.specs, FUN = readFwf)
    loaded.data <- lapply(extended.specs, FUN = prepareRetval)

    return(loaded.data)
}
