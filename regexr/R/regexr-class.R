#' Prints a regexr Object
#' 
#' Prints a \code{regexr} object
#' 
#' @param x The \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @method print regexr
print.regexr <- function(x, ...){
    print(as.character(x))
}

#' Unglue regexr Object
#' 
#' \code{regexr} - unglue \code{regexr} object.
#' 
#' regexr Method for unglue
#' @param x The \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @method unglue regexr
unglue.regexr <- function(x, ...){
    out <- attributes(x)[["subs"]]
    class(out) <- unique(c("unglued", class(out)))
    out
}

#' Prints an unglued object
#' 
#' Prints an unglued object.
#' 
#' @param x The unglued object.
#' @param \ldots Ignored.
#' @export
#' @method print unglued
print.unglued <- function(x, ...){
    class(x) <- "list"
    out <- invisible(lapply(x, function(x) {class(x) <- "character";x}))
    print(out)
}



#' Extract Comments From regexr Object
#' 
#' \code{regexr} - Extract comments from \code{regexr} object.
#' 
#' regexr Method for comments
#' @param x The \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @method comments regexr
comments.regexr <- function(x, ...){
    attributes(x)[["comments"]]
}

#' Set Comments For regexr
#' 
#' \code{regexr} - Set comments for \code{regexr} object.
#' 
#' regexr Method for comments<-
#' @param x The \code{regexr} object.
#' @param value The comment(s) to assign.
#' @param \ldots Ignored.
#' @export
#' @method comments<- regexr
`comments<-.regexr` <- function(x, value){
    attributes(x)[["comments"]] <- value
    len <- length(attributes(x)[["comments"]])
    dif <- diff(c(length(attributes(x)[["subs"]]), len))
    if (dif > 0){  
        null <- structure(list(NULL), .Names = "")
        attributes(x)[["subs"]] <- unlist(list(attributes(x)[["subs"]], 
            rep(null, dif)), recursive=FALSE)
    }    
    x
}

#' Get Sub-expressions From a regexr Object
#' 
#' \code{subs} - Get the sub-expressions from a \code{regexr} object.
#' 
#' regexr Method for subs
#' @param x The \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @method subs regexr
subs.regexr <- function(x, ...){
    attributes(x)[["subs"]]
}

#' Set Regex Sub-expressions From a regexr Object
#' 
#' \code{subs<-} - Set the sub-expressions of a \code{regexr} object.
#' 
#' regexr Method for subs<-
#' @param x The \code{regexr} object.
#' @param value The comment(s) to assign.
#' @param \ldots Ignored.
#' @export
#' @method subs<- regexr
`subs<-.regexr` <- function(x, value){
    attributes(x)[["subs"]] <- value
    len <- length(attributes(x)[["subs"]])
    dif <- diff(c(length(attributes(x)[["comments"]]), len))
    if (dif > 0){  
        null <- structure(list(NULL), .Names = "")
        attributes(x)[["comments"]] <- unlist(list(attributes(x)[["comments"]], 
            rep(null, dif)), recursive=FALSE)
    }
    x[[1]] <- paste(unlist(attributes(x)[["subs"]]), collapse="")
    x
}




#' Get Names of Sub-Expressions of a regexr Object
#' 
#' Get names of a \code{regexr} object.
#' 
#' @param x The \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @method names regexr
names.regexr <- function(x, ...){

    names(attributes(x)[["subs"]])

}

#' Set Names of a Sub-expressions of a regexr Object
#' 
#' Set names of a \code{regexr} object's sub-expressions.
#' 
#' @param x The \code{regexr} object.
#' @param value The comment(s) to assign.
#' @param \ldots Ignored.
#' @export
#' @method names<- regexr
`names<-.regexr` <- function(x, value){

    rnull <- is.null(names(attributes(x)[["subs"]]))
    cnull <- is.null(names(attributes(x)[["comments"]]))

    names(attributes(x)[["subs"]]) <- value
    if (rnull) {
        names(attributes(x)[["subs"]])[is.na(names(attributes(x)[["subs"]]))] <- ""
    }
    names(attributes(x)[["comments"]]) <- value
    if (cnull) {
        names(attributes(x)[["comments"]])[is.na(names(attributes(x)[["comments"]]))] <- ""
    }
    x

}

#' Summarize a regexr Object
#' 
#' Summarize a \code{regexr} object.
#' 
#' @param object The \code{regexr} object 
#' @param \ldots Ignored.
#' @method summary regexr
#' @export
summary.regexr <- function(object, ...){

    if (length(attributes(object)[["comments"]]) != 
        length(attributes(object)[["subs"]])) {
        warning("Mismatch in number of subs and comments; items recycled\n",
            "Consider using `comments` and/or `subs` to update the regexr object")
    }
    out <- suppressWarnings(Map(function(x, y) list(comment = x, subs=y),  
        attributes(object)[["comments"]],
        attributes(object)[["subs"]]
    ))
    class(out) <- "summary_regexr"
    attributes(out)[["subs"]] <- as.character(object)
    out
}


#' Prints a summary_regexr object
#' 
#' Prints a summary_regexr object.
#' 
#' @param x The summary_regexr object.
#' @param \ldots Ignored.
#' @export
#' @method print summary_regexr
print.summary_regexr <- function(x, ...){
   
    class(x) <- "list"

    reg <- attributes(x)[["subs"]]
    cat("\n", reg, "\n", 
        paste(rep("=", nchar(reg)), collapse=""), "\n"
    )

    x <- namer(x)

    for (i in seq_along(x)) {
        element <- sprintf("SUB-EXPR %s: ", i)
        len <- nchar(element)
        message(element,  x[[i]][["subs"]])
        message(sprintf("NAME%s: ", paste(rep(" ", len - 6), collapse="")),  names(x)[i])        
        message(sprintf("COMMENT%s: ", paste(rep(" ", len - 9), collapse="")), 
            sprintf("\"%s\"", x[[i]][["comment"]]), "\n")
    }
}

#' Test Regular Expression Validity
#' 
#' Test regular expression validity of a \code{regexr} object.
#' 
#' test Method for subs<-
#' @param x The \code{regexr} object.
#' @param quiet logical.  Should \code{test} print warnings about the 
#' concatenated expression and individual sub-expressions?
#' @param \ldots Ignored.
#' @export
#' @method test regexr
test.regexr <- function(x, quiet = FALSE, ...){

    out1 <- is.regex(x)
    if (!isTRUE(quiet) && !out1){
        warning("The concatenated regex is not valid\n\n", as.character(x), "\n")
    }
    out2 <- sapply(subs(x), is.regex)
    if (!isTRUE(quiet) && any(!out2)){
        warning("The following regex sub-expressions are not valid in isolation:\n\n",
            paste(paste0("(", seq_len(sum(!out2)), ") ", 
                as.character(unlist(subs(x)))[!out2]), collapse="\n")
        )
    }
    list(regex = out1, subexpressions = out2)
}


