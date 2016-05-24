## Data tables are simplified data frames without coercion.
## All variables in the data table must be vectors of the same length.
## Note that in current versions of R, data frame conversions should be
## confined to creation via data.frame(), but not happen in subsequent
## computations:
##
##   data.frame() converts each of its arguments to a data frame by
##   calling as.data.frame(optional=TRUE).  As that is a generic
##   function, methods can be written to change the behaviour of 
##   arguments according to their classes: R comes with many such
##   methods.  Character variables passed to data.frame() are converted 
##   to factor columns unless protected by I().  If a list or data
##   frame or matrix is passed to data.frame() it is as if each
##   component or column had been passed as a separate argument (except
##   for matrices of class 'model.matrix' and those protected by I()).
##
## Arguably, it would be nicer if data.frame() simply put its arguments
## into a list of variables, which is what data_table() does.  But then
## the low level representation dictates what happens.  E.g., one cannot
## put POSIXlt objects into a data table because these are implemented
## as lists of the respective date/time components.  Oh well.

## We currently implement data tables as data frame like objects.
## The data_table() creator puts all given arguments into the table; row
## names need to be set lateron.  All arguments must be named.
## Package relations also has an as.data.frame.relation() method taking
## an optional row.names argument.  Hence, we have an internal workhorse
## .make_data_frame_from_list() which validates row names against the
## data, but does not check argument lengths, the equality of which
## needs to be ensured by the context or the caller.

## <FIXME>
## How should names be handled?  Via deparsing as for data.frame()?
## Should unique names be ensured?
## </FIXME>

data_table <-
function(...)
{
    args <- list(...)
    classes <- c("data_table", "data.frame")
    if(!length(args))
        return(`class<-`(data.frame(), classes))
    if(is.null(nms <- names(args)) || !all(nzchar(nms)))
        stop("All arguments must be named.")
    ## Check lengths.
    if(length(table(sapply(args, length))) > 1L)
        stop("All arguments must have the same length.")
    `class<-`(.make_data_frame_from_list(args), classes)
}

print.data_table <-
function(x, ...)
{
    if(!length(x)) {
        writeLines(gettextf("An empty data table."))
        return(invisible(x))
    }
    ## Should this try to use a row.names = FALSE default?
    NextMethod("print")
}
