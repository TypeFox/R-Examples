#' Replace null values
#' 
#' @param x Value to be tested.
#' @param replacement Value to use if \code{x} is \code{NULL}.
#' @return \code{x} or \code{replacement} in case \code{x} is \code{NULL}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
ifnull <- function(x, replacement=NA) if(is.null(x)) replacement else x

#' Display contents of an evironment, data.frame or list as a summary table
#'
#' Color coded according to class and dimensions of contents. See
#' \code{\link[xtermStyle]{style}} for details.
#'
#' @param envir Environment, data frame or list to be displayed. Optional,
#'   default: globalenv()
#' @param pattern Regexp filtering of objects. Only objects matching the pattern
#'   are displayed. Optional, default: show all objects.
#' @param all.names Whether to show hidden objects.
#' @param exclude A list of objects not to be displayed. To set a default exclusion
#'   mask use the \code{whos.set.mask} function. If \code{whos.set.mask} is
#'   called without a list of object names all objects currently in globalenv()
#'   are hidden. This is useful for example if you have a lot of stuff in the
#'   workspace that you aren't currently interested in but is needed to make
#'   your code run.
#' @return Nothing
#' @examples
#' whos()
#' data(USArrests)
#' whos(USArrests)
#' 
#' data(iris)
#' whos()
#' whos.all()
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @import data.table
#' @import xtermStyle
#' @seealso whos.options, browse
#' @export
whos <- function(envir=parent.frame(), pattern=".", all.names, exclude=getOpt("exclude")){
    # Interpret the `envir` argument if not already an environment
    if(is.function(envir)){
        envir <- environment(envir)
    } else {
        envir <- switch(class(envir)[1],
            `character` = {
                pattern <- envir
                parent.frame()
            },
            `integer` = as.environment(envir),
            `numeric` = {
                if(envir != as.integer(envir)) stop("Invalid `envir` argument")
                as.environment(envir)
            },
            envir
        )
    }
    
    # Get names of objects in environment matching pattern
    accessors <- if(is.environment(envir)){
        if(missing(all.names)){
            l <- ls(envir=envir, all.names=FALSE)
            if(length(l) == 0)
                l <- ls(all.names=TRUE, envir=envir)
        } else {
            l <- ls(all.names=all.names, envir=envir)
        }
        structure(l, names=l)
    } else if(isS4(envir)){
        structure(slotNames(envir), names=slotNames(envir))
    } else {
        # lists, data.frames, data.tables etc
        if(!is.null(names(envir))){
            structure(seq_len(length(envir)), names=names(envir))
        } else {
            1:length(envir)
        }
    }
    if(!is.null(names(accessors))){
        accessors <- accessors[!names(accessors) %in% exclude &
                               grepl(pattern, names(accessors))]
    }

    if(length(accessors) == 0){
        NULL
    } else {
        obj.sapply <- if(isS4(envir)){
            function(fun, ...) sapply(accessors, function(x) ifnull(fun(slot(envir, x))), ...)
        } else if(is.character(accessors)){
            # Unordered set of objects, as found in an environment
            function(fun, ...) sapply(accessors, ..., FUN = function(x){
                ifnull(tryCatch(fun(get(x, envir)), error=function(...) NULL))
            })
        } else {
            # Ordered set of objects, as found in a data frame
            # Use the index numbers rather than names in case of duplicates
            function(fun, ...) sapply(envir, function(x, ...) ifnull(fun(x, ...)), ...)
        }

        # Make an object/property matrix (objects as rows, properties as columns)
        env.cols <- lapply(getOpt("columns")$envir, function(f) f(envir, accessors))
        structure(do.call(data.table, c(
            list(
                style = obj.sapply(style.auto),
                name = if(!is.null(names(accessors))) names(accessors) else NA
            ),
            env.cols[!sapply(env.cols, is.null)],
            lapply(getOpt("columns")$object, obj.sapply)
        )), class=c("whos", "data.table", "data.frame"))
    }
}

#' @param x \code{\link{whos}} object.
#' @noRd
#' @export
print.whos <- function(x, ...){
    # Remove columns without content
    # Always keep the style column, even if no style is used, since it has a special function
    x <- x[, c(TRUE, !sapply(x, function(x) all(x %in% c(FALSE, NA)))[-1]), with=FALSE]

    # Calculate summaries
    xnames <- names(x)
    xsum <- mapply(function(field, values){
        fun <- ifnull(getOpt("summary")[[field]], function(...) NA)
        val2str <- ifnull(getOpt("print")[[field]],
                          function(x) if(is.na(x)) "" else x)
        val2str(fun(values))
    }, names(x), x)

    # Convert all to characters
    x <- matrix(ncol=ncol(x), mapply(function(field, values){
        val2str <- getOpt("print")[[field]]
        if(is.null(val2str)){
            if(is.logical(values)){
                ifelse(is.na(values), "", ifelse(values, sprintf("[%s]", field), ""))
            } else ifelse(is.na(values), "", as.character(values))
        } else val2str(values)
    }, names(x), x))

    nc <- pmax(nchar(xnames), apply(x, 2, function(x) max(nchar(x))))[-1]
    align <- ifelse(getOpt("align")[xnames[-1]] %in% "right", "", "-")
    fmt <- paste("%s ", paste(sprintf("%%%s%is", align, nc), collapse="  "), " ", style.clear(), "\n", sep="")
    tryCatch({
        index.nchar <- ceiling(log10(1+nrow(x)))

        # Field names
        cat(sprintf(sprintf("%%%is", index.nchar+2), ""),
            do.call(sprintf, as.list(c(fmt, "", xnames[-1]))), sep="")

        # Content
        cat(paste(
            sprintf(sprintf("%%%ii: ", index.nchar), seq_len(nrow(x))),
            apply(x, 1, function(x) do.call(sprintf, as.list(c(fmt, x)))),
            sep="", collapse=""))
        
        # Summary
        if(any(xsum != "") && nrow(x) > 1){
            cat("Summary:\n",
                sprintf(sprintf("%%%is", index.nchar+2), ""),
                do.call(sprintf, as.list(c(fmt, xsum))), sep="")
        }
    }, interrupt=cat(style.clear()))
}
#' Subset a whos object
#'
#' This is exactly the same as the subsetting function of data.table, but
#' preserves the class, to make it find the right print function.
#'
#' @param x \code{\link{whos}} object.
#' @param ... Sent to 
#' @noRd
#' @export
`[.whos` <- function(x, ...){
    as.whos(as.data.table(x)[...])
}

#' Set default behavior of the whos function
#'
#' @param exclude Objects to exclude from view. Can be a character vector of
#'   names, or an environment, but not any regular expressions so far.
#' @param report.S4.size Calculating the size of S4 objects with
#'   \code{\link{object.size}} can take an annoyingly long time (seconds), set
#'   this option to \code{FALSE} to skip it and get quicker execution.
#' @return Nothing. The values are stored as global options.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
whos.options <- function(exclude, report.S4.size){
    if(!missing(exclude)){
        options(whos.exclude = switch(class(exclude),
            `character` = exclude,
            `environment` = ls(exclude),
            `integer` = ls(as.environment(exclude))
        ))
    }
    if(!missing(report.S4.size))
        options(whos.report.S4.size = report.S4.size)
}

#' @param x A character vector of object names to exclude or include.
#' @param pattern \link[=regex]{Regular expression pattern} to match object
#'   names against, e.g. \code{pattern="^my\\..*"} will exclude or include
#'   \code{"my.vector"} and \code{"my.matrix"} but not \code{"mysql.con"}.
#' @param envir Environment to search in.
#' @rdname whos.options
#' @export
whos.exclude <- function(x=NULL, pattern, envir=parent.frame()){
    if(!missing(pattern)) x <- union(x, ls(envir=envir, pattern=pattern))
    p <- getOption("dataview")
    p$exclude <- union(p$exclude, x)
    options(dataview = p)
}
#' @rdname whos.options
#' @export
whos.include <- function(x=NULL, pattern, envir=parent.frame()){
    if(!missing(pattern)) x <- union(x, ls(envir=envir, pattern=pattern))
    p <- getOption("dataview")
    p$exclude <- setdiff(p$exclude, x)
    options(dataview = p)
}

#' Shortcut for calling whos without exclusion.
#'
#' @param ... Parameters sent to \code{\link{whos}}.
#' @return Nothing
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @rdname whos
#' @export
whos.all <- function(...){
    whos(..., exclude=NULL)
}

#' @param x \code{\link{whos}} object.
#' @param keep.rownames Ignored, kept for S3 consistency.
#' @noRd
#' @export
as.data.table.whos <- function(x, keep.rownames=FALSE){
    class(x) <- setdiff(class(x), "whos")
    x
}

#' @param x \code{\link{whos}} object.
#' @param optional Ignored, kept for S3 consistency.
#' @param row.names Ignored, kept for S3 consistency.
#' @param ... Ignored, kept for S3 consistency.
#' @noRd
#' @export
as.data.frame.whos <- function(x, row.names=NULL, optional=FALSE, ...){
    class(x) <- setdiff(class(x), c("whos", "data.table"))
    x
}

#' Convert objects to whos
#'
#' @param x Object of type \code{\link{data.table}} or \code{\link{data.frame}}.
#' @examples
#' an.object <- "Containing all my stuff"
#' w <- as.data.frame(whos())
#' as.whos(w)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
as.whos <- function(x){
    UseMethod("as.whos")
}
#' @rdname as.whos
#' @export
as.whos.data.table <- function(x){
    class(x) <- c("whos", class(x))
    x
}
#' @rdname as.whos
#' @export
as.whos.data.frame <- function(x){
    as.whos(as.data.table(x))
}
