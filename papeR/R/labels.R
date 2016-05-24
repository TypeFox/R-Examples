################################################################################
##  Author: Benjamin Hofner, benjamin.hofner@fau.de

################################################################################
# Extract labels from data sets
labels.data.frame <- function(object, which = NULL, abbreviate = FALSE, ...) {

    ## if no labels were properly set use alternative methods to specify labels:
    if (!is.ldf(object)) {
        ## if no labels specified temporarily set names as labels
        if (is.null(attr(object, "variable.labels"))) {
            labels(object) <- colnames(object)
        } else {
            ## clean labels
            object <- CLEAN_LABELS(object)
            ## set these labels temporarily as elementwise labels
            labels(object) <- attr(object, "variable.labels")
        }
    }

    ## which labels should be extracted?
    which <- check_which(which, object, "extract")

    ## now extract labels
    RET <- sapply(as.data.frame(object)[which], get_labels)
    ## fix non-existing labels
    if (is.list(RET) && any(idx_null <- sapply(RET, is.null))) {
        nms <- colnames(object)
        if (is.character(which))
            names(nms) <- nms
        RET[idx_null] <- nms[which][idx_null]
        RET <- unlist(RET)
    }

    ## should labels be abbreviated?
    if (abbreviate) {
        nms <- names(RET)
        RET <- abbreviate(RET, ...)
        names(RET) <- nms
    }
    return(RET)
}

################################################################################
# Extract labels from labeled variables
labels.lv <- function(object, abbreviate = FALSE, ...) {
    RET <- get_labels(object)
    ## should labels be abbreviated?
    if (abbreviate) {
        nms <- names(RET)
        RET <- abbreviate(RET, ...)
        names(RET) <- nms
    }
    RET
}

################################################################################
# Sets labels
"labels<-" <- function(data, which = NULL, value){

    which <- check_which(which, data, "define")

    if (!is.null(value)) {
        if (length(which) != length(value))
            stop("One must supply a label for each _selected_ column of the data set.")
        if (is.character(which))
            names(value) <- which
    }

    for (i in which) {
        attr(data[[i]], "variable.label") <- value[[i]]
        class(data[[i]]) <- c("lv", class(data[[i]]))
    }

    ## remove attribute of data set if it exists
    if (!is.null(attr(data, "variable.labels")))
        attr(data, "variable.labels") <- NULL

    class(data) <- c("ldf", class(data))
    return(data)
}

"labels[<-" <- function(data, i, value)
    labels(data, which = i) <- value


CLEAN_LABELS <- function(data) {
    ## drop spare labels
    spare <- !(names(attr(data, "variable.labels")) %in% names(data))
    if (any(spare)) {
        message("Note: Variables have been removed or label names and ",
                "column names don't match. ",
                "Corresponding variable labels are removed.")
        attr(data, "variable.labels") <-  attr(data, "variable.labels")[!spare]
    }
    ## add missing labels
    missing <- !(names(data) %in% names(attr(data, "variable.labels")))
    if (any(missing)) {
        tmp <- names(data)[missing]
        names(tmp) <- names(data)[missing]
        attr(data, "variable.labels") <- c(attr(data, "variable.labels"),
                                           tmp)
    }
    ## re-order
    attr(data, "variable.labels") <- attr(data, "variable.labels")[names(data)]
    ## return altered data set
    return(data)
}


## define coercion function
as.ldf <- function(object, ...)
    UseMethod("as.ldf")

as.ldf.data.frame <- function(object, ...) {
    labels(object) <- labels(object)
    object
}

convert.labels <- function(object)
    as.ldf.data.frame(object)

is.ldf <- function(object)
    !all(sapply(lapply(object, get_labels), is.null))
