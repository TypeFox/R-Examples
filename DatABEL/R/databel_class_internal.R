#
# databel internal util
#

databel_check <- function(x, reconnect = TRUE, stop_on_error = TRUE,
                          quiet = FALSE)
{
    #	print("dataprovcbel_check started");

    if (class(x) != "databel") {
        msg <- "databel_check: object is not of class 'databel'"
        stop(msg)
        #if (stop_on_error) {stop(msg);} else {warning(msg);return(FALSE)}
    }

    if (class(x@data) != "externalptr") {
        msg <- "databel_check: data is not of class 'externalptr'";
        stop(msg);
        #if (stop_on_error) {stop(msg);} else {warning(msg);return(FALSE)}
    }

    if (externalptr_is_null(x@data)) {
        if (!reconnect) {
            msg <- "databel_check: object is not connected (will not work for writing, access is slower); use 'connect(object)'"
            if (stop_on_error) {
                stop(msg);
            } else {
                if (!quiet) warning(msg);
                return(FALSE)
            }
        } else {
            result <- try(
                eval.parent(substitute(connect(x)))
                )
            if (class(result) != "try-error") {
                return(TRUE)
            } else {
                msg <- "databel_check: can not connect object"
                if (stop_on_error) {stop(msg);} else {warning(msg);return(FALSE)}
            }
        }
    }
    #	print("databel_check finished");
    return(TRUE);
};


externalptr_is_null <- function(x) {
    if (!is(x, "externalptr")) stop("x is not 'externalptr'")
    return(.Call("externalptr_is_null", x))
}


uninames <- function(filtredmatrixptr)
{
    out <- list()
    out$unique.names <- out$unique.colnames <- out$unique.rownames <- FALSE
    colnames <- .Call("get_all_varnames_R", filtredmatrixptr)
    rownames <- .Call("get_all_obsnames_R", filtredmatrixptr)

    if (!anyDuplicated(colnames)) {
        out$unique.colnames <- TRUE
    } else {
        warning("uninames: some column names are not unique; use set_dimnames/get_dimnames for non-unique row/col names")
    }

    if (!anyDuplicated(rownames)) {
        out$unique.rownames <- TRUE
    } else {
        warning("uninames: some row names are not unique; use set_dimnames/get_dimnames for non-unique row/col names")
    }

    if (all(colnames == c(1:length(colnames)))) out$unique.colnames <- FALSE
    if (all(rownames == c(1:length(rownames)))) out$unique.rownames <- FALSE

    out$unique.names <- (out$unique.colnames && out$unique.rownames)

    return(out)
}
