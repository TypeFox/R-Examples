#' Retrieves package options
#' 
#' The dataview package contains a number of option for tailoring its
#' behaviour. These options are stored as a named list in a single global
#' varible named "dataview". To overwrite an option fetch the default values
#' with \code{default.options}, modify the returned results and set it back as
#' in the example below.
#'
#' Below is a description of the available options,
#' but to understand how they work it is probably easier to directly study the
#' return value of \code{default.options}.
#' \describe{
#'     \item{\code{align}}{Column alignment, left or right.}
#'     \item{\code{columns}}{A named list with columns to use in \code{\link{whos}}.
#'         The \code{object} element should be a named list of functions to be
#'         applied on each object that is to be queried. 
#'         The \code{envir} element should be a named list of functions to be
#'         applied on the object or environment on which \code{\link{whos}} is
#'         called together with its "accessors", which are object names
#'         (for environments) or named indices (for everything else).
#'         All functions should return a vector of values that can be used as
#'         a \code{\link{data.table}} column.}
#'     \item{\code{print}}{A named list of custom print functions.
#'         By default each column of a \code{\link{whos}} object
#'         is printed in a similar way as \code{\link{data.table}} or
#'         \code{\link{data.frame}}. If you wish to override this behaviour for
#'         a given column, please supply a named print function here.
#'         The function will be given the column as returned by the corresponding
#'         column function above, and should produce a character vector.}
#'     \item{\code{summary}}{A list of named summary functions.
#'         These are fed a column and should return a single value.}
#' }
#' @param x Option to retrieve.
#' @examples
#' # This shows how to modify a column presented by whos.
#' # The new function only reports the size of non-S4 objects
#' # to improve execution time.
#' opt <- default.options()
#' opt$columns$bytes <- function(x) if(isS4(x)) NA else object.size(x)
#' options(dataview = opt)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
getOpt <- function(x){
    ifnull(getOption("dataview")[[x]],
           default.options()[[x]])
}
#' @rdname getOpt
#' @export
default.options <- function(){
    list(
        align = c(bytes="right"),
        columns = list(
            envir = list(
                Key = function(env, nam){
                    if(!is.data.table(env) || !haskey(env)) return(NULL)
                    names(nam) %in% key(env)
                },
                argument = function(env, nam){
                    if(!is.environment(env) || identical(env, .GlobalEnv))
                        return(NULL)
                    stack.level <- sapply(sys.frames(), identical, env)
                    if(!any(stack.level)) return(NULL)
                    f <- formals(sys.function(which(stack.level)))
                    if(is.null(f)) return(NULL)
                    f.modified <- sapply(names(f), function(x){
                        tryCatch(!identical(f[[x]], get(x, env)),
                                 error=function(err) NA)
                    })
                    f.missing <- sapply(f, identical, quote(expr=))
                    # Missing arguments are represented by the "empty symbol",
                    # also produced by quote(expr=).
                    # Read more at http://stackoverflow.com/a/27824791/840460
                    factor(ifelse(!nam %in% names(f), NA,
                           ifelse(!is.na(f.modified[nam]) & f.modified[nam], "custom",
                           ifelse(f.missing[nam], "missing", "default"))),
                           c("default", "custom", "missing"))
                }
            ),
            object = list(
                class = function(x){
                    cls <- class(x)
                    paste0(cls[1],
                           if(length(cls) > 1) sprintf(" (+%i)", length(cls)-1)
                           else "", sep="")
                },
                S4 = isS4,
                dim = function(x){
                    if(is.function(x) || is.null(x)){
                        NA
                    } else {
                        x <- list(length = length(x), dim = dim(x))
                        if(is.null(x$dim)){
                            as.character(x$length)
                        } else {
                            paste(x$dim, collapse="x")
                        }
                    }
                },
                bytes = object.size,
                comment = function(x) !is.null(comment(x))
            )
        ),
        print = list(
            bytes = function(x){
                unit <- c("B  ", "B  ", "KiB", "MiB", "GiB", "TiB", "EiB")[
                    1+sapply(x, function(b) sum(b > 1024^(0:4)))]
                size <- ifelse(x == 0, 0, 2^(log2(x) %% 10))
                ifelse(is.na(x), "", sprintf("%.4g %s", size, unit))
            },
            comment = function(x) ifelse(x, "+!", "")
        ),
        summary = list(bytes = sum)
    )
}
