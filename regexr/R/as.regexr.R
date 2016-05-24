#' Generic Method to Coerce to regexr
#' 
#' Coerce an object to \code{regexr} class.
#' 
#' @param x An object to coerce to a \code{regexr} object.
#' @param names logical.  Should names be included in the \code{construct} 
#' script?
#' @param comments  logical.  Should comments be included in the \code{construct} 
#' script?
#' @param names.above logical.  Should ames be included above the regex
#' in the \code{construct} script?  If \code{FALSE} names are placed in front of 
#' the sub-expressions.
#' @param comments.below  logical.  Should comments be included below the 
#' sub-expressions in the \code{construct} script?  If \code{FALSE} comments 
#' are placed behind the sub-expressions.
#' @param \ldots Other arguments passed to \code{as.regexr} methods.
#' @return Returns a dual \code{regexr} and \code{reverse_construct} object.
#' @export
#' @note \code{as.regexr.character} not currently supported as it utilized 
#' 'http://rick.measham.id.au/paste/explain', however this website appears not 
#' to be operational.  This functionality may return if the website is again 
#' available or if a suitable replacement regex parser is found.  
#' 
#' If you have a suggested web based parser please make a suggestion via:
#' \url{https://github.com/trinker/regexr/issues}.
## @note \code{as.regexr.character} utilizes \url{http://rick.measham.id.au/paste/explain}
## to break the regular expression into sub-expressions.
#' @examples
#' \dontrun{
#' ## NOTE: These examples will likely fail unless 
#' ## http://rick.measham.id.au/paste/explain becomes available again.
#' library("qdapRegex")
#' (myregex <- grab("@@rm_time2"))
#' out <- as.regexr(myregex)
#' 
#' out
#' summary(out)
#' comments(out)
#' subs(out)
#' test(out)
#' get_construct(out)
#' 
#' ## On Windows copy to clipboard
#' get_construct(out, file="clipboard")
#' 
#' ## No names & comments behind sub-expressions
#' myregex2 <- "(\\s*[a-z]+)([^)]+\\))"
#' get_construct(as.regexr(myregex2, names=FALSE))
#' get_construct(as.regexr(myregex2, names=FALSE, names.above = TRUE, 
#'     comments.below = TRUE))
#' }
as.regexr <- function(x, names = TRUE, comments = TRUE, names.above = FALSE, 
    comments.below = FALSE, ...){
    UseMethod("as.regexr")
}


#' Coerce character to regexr
#' 
#' Convert a regular expression to a commented \code{regexr} object.
#' 
#' character Method for as.regexr
#' @param x The \code{character} object.
#' @param names logical.  Should names be included in the \code{construct} 
#' script?
#' @param comments  logical.  Should comments be included in the \code{construct} 
#' script?
#' @param names.above logical.  Should ames be included above the sub-expressions
#' in the \code{construct} script?  If \code{FALSE} names are placed in front of 
#' the sub-expressions.
#' @param comments.below logical.  Should comments be included below the 
#' sub-expressions in the \code{construct} script?  If \code{FALSE} comments are 
#' placed behind the sub-expressions.
#' @param \ldots Ignored.
#' @export
#' @method as.regexr character
as.regexr.character <- function(x, names = TRUE, comments = TRUE, 
    names.above = FALSE, comments.below = FALSE, ...){

    out <- try(regex_break_down(x))

    if (inherits(out, "try-error")) { # added 8-16-15 when http://rick.measham.id.au/paste went down
        message("The `http://rick.measham.id.au/paste/explain` regex conversion no longer", 
            "\nworks as the website currently does not work.\n\n",
            "The regex conversion functionality will return if the website becomes",
            "\noperational again or if a suitable substitute can be found."
        )
        return(NULL)
    }
    
    loc <- gregexpr("EXPLANATION", out[1])
    out <- out[-c(1:2)]

    breaks <- grepl("^-{10,}$", out)
    inds <- !breaks

    out <- split(out[inds], cumsum(breaks)[inds])
    names(out) <- as.numeric(names(out)) + 1

    pieces <- lapply(out, function(x){
        y1 <- gsub("\\s+$", "", substring(x, 1, loc))
        y1[-1] <- gsub("^\\s+", "", y1[-1])
        y2 <- gsub("^\\s+|\\s+$", "", substring(x, loc))
        lets <- c("n", "r", "t", "f", "a")
        for (i in seq_len(length(lets))){
            y2 <- gsub(paste0("\\\\", lets[i]), paste0("\\", lets[i]), y2, fixed=TRUE)
        }
        y1 <- paste(y1, collapse="")
        y2 <- paste(y2, collapse=" ")

        list(regex = y1, comment = y2)
    })

    pieces4regexr <- lapply(pieces, function(x){
        x[[1]] <- gsub("\\\\", "\\", gsub("^\\s+", "", x[[1]]), fixed=TRUE)
        x
    })

    out <- x

    class(out) <- c("regexr", "reverse_construct", class(out))
    attributes(out)[["subs"]] <- stats::setNames(sapply(pieces4regexr, "[", 1), 
        names(pieces4regexr))
    attributes(out)[["comments"]] <- stats::setNames(sapply(pieces4regexr, "[", 2), 
        names(pieces4regexr))

    if (!comments.below) {
        max.nchar.regex <- max(sapply(pieces, function(x) nchar(x[[1]])))
    }

    pieces4construct <- sapply(seq_along(pieces), function(i){
        x <- pieces[[i]]
        x[[1]] <- gsub("\"", "\\\\\"", x[[1]])
        x[[2]] <- gsub("\"", "\\\\\"", x[[2]])
        x[[1]] <- gsub(" ", "  ", x[[1]])
        indent <- (nchar(x[[1]]) - nchar(gsub("^\\s+", "", x[[1]]))) 
        x[[1]] <- paste0("    ", x[[1]])

        if (isTRUE(names)) {
            if (isTRUE(names.above)) {
                thenames <- paste0(paste(rep(" ", indent), collapse=""), "`", 
                    names(pieces)[i], "` = \n")
            } else {
                thenames <- paste0(paste(rep(" ", indent), collapse=""), "`", 
                    names(pieces)[i], "` = ")
            }
        } else {
            thenames <- ""
        }
        
        if (isTRUE(names)) {
            if (isTRUE(names.above)) {
                theregexes <- gsub("(^\\s+)", "\\1\"", x[[1]])
            } else {
                theregexes <- gsub("^\\s+", "\\1\"", x[[1]])
            }
        } else {
            theregexes <- gsub("^\\s{4}", "", gsub("(^\\s+)", "\\1\"", x[[1]]))
        }

        if (isTRUE(comments)) {
            if (isTRUE(comments.below)) {
                thecomments <- paste0("\"\n", paste(rep(" ", indent + 8), 
                    collapse=""), "%:)%\"", x[[2]], "\"")
            } else {
                thecomments <- paste0("\"", 
                    paste(rep(" ", indent + (max.nchar.regex - nchar(theregexes)) + 10), 
                    collapse=""), "%:)%  \"", x[[2]], "\"")
            }
        } else {
            thecomments <- ""
        }
        paste0(thenames, theregexes, thecomments)
    })

    reverse_construct <- paste0("construct(\n", 
        paste(pieces4construct, collapse=",\n"), "\n)\n")
    class(reverse_construct) <- c("reverse_construct", class(reverse_construct))
    attributes(out)[["reverse_construct"]] <- reverse_construct

    out
}


#' Prints a reverse_construct object
#' 
#' Prints a reverse_construct object.
#' 
#' @param x A \code{reverse_construct} object.
#' @param file A connection, or a character string naming the file to print to. 
#' If "" (the default), \code{\link[base]{cat}} prints to the console unless 
#' redirected by \code{\link[base]{sink}}.  Windows users may use 
#' \code{file = "clipboard"} to copy the content to the clipboard.  
#' @param \ldots Other arguments passed to \code{\link[base]{cat}}.
#' @export
#' @method print reverse_construct
print.reverse_construct <- function(x, file = "", ...){
    cat(x, file = file, ...)
}



#' Coerce default to regexr
#' 
#' Convert a regular expression to a commented \code{regexr} object.
#' 
#' default Method for as.regexr
#' @param x The object to be corced to \code{regexr}.
#' @param names logical.  Should names be included in the \code{construct} 
#' script?
#' @param comments  logical.  Should comments be included in the \code{construct} 
#' script?
#' @param names.above logical.  Should ames be included above the sub-expressions
#' in the \code{construct} script?  If \code{FALSE} names are placed in front of 
#' the sub-expressions.
#' @param comments.below  logical.  Should comments be included below the 
#' sub-expressions in the \code{construct} script?  If \code{FALSE} comments are 
#' placed behind the sub-expressions.
#' @param \ldots Ignored.
#' @export
#' @method as.regexr default
as.regexr.default <- as.regexr.character


#' Extract Script from \code{reverse_construct} to \code{construct} a 
#' \code{regexr} Object.
#' 
#' Pulls the \code{reverse_construct} attribute from a \code{reverse_construct}.  
#' This script can be assigned to an object and run in the console to create a 
#' comments, named \code{regexr} object.
#' 
#' @param x A \code{reverse_construct} object.
#' @param file A connection, or a character string naming the file to print to. 
#' If "" (the default), \code{\link[base]{cat}} prints to the console unless 
#' redirected by \code{\link[base]{sink}}.  Windows users may use 
#' \code{file = "clipboard"} to copy the content to the clipboard.
#' @param \ldots Other arguments passed to \code{\link[regexr]{print.reverse_construct}}.
#' @return Returns an auto-commented script used to \code{construct} a 
#' \code{regexr} object.
#' @export
#' @examples
#' \dontrun{
#' library("qdapRegex")
#' (myregex <- grab("@@rm_time2"))
#' out <- as.regexr(myregex)
#' 
#' out
#' summary(out)
#' comments(out)
#' subs(out)
#' test(out)
#' get_construct(out)
#' 
#' ## On Windows copy to clipboard
#' get_construct(out, file="clipboard")
#' }
get_construct <- function(x, file = "", ...){
    UseMethod("get_construct")
}

#' Extract Script from \code{reverse_construct} to \code{construct} a 
#' \code{regexr} Object.
#' 
#' Pulls the \code{reverse_construct} attribute from a \code{reverse_construct}.  
#' This script can be assigned to an object and run in the console to create a 
#' comments, named \code{regexr} object.
#' 
#' reverse_construct Method for get_construct
#' @param x A \code{reverse_construct} object.
#' @param file A connection, or a character string naming the file to print to. 
#' If "" (the default), \code{\link[base]{cat}} prints to the console unless 
#' redirected by \code{\link[base]{sink}}.  Windows users may use 
#' \code{file = "clipboard"} to copy the content to the clipboard. To print
#' as \code{character} use \code{file = NULL}.
#' @param \ldots Other arguments passed to \code{\link[regexr]{print.reverse_construct}}.
#' @return Returns an auto-commented script used to \code{construct} a 
#' \code{regexr} object.
#' @export
#' @method get_construct reverse_construct
get_construct.reverse_construct <- function(x, file = "", ...){

    out <- attributes(x)[["reverse_construct"]]
    if (!is.null(file)) {
        print(out, file = file, ...)
    } else {   
        as.character(out)
    }    

}




regex_break_down <- function(pattern){
   
    URL2 <- paste0("http://rick.measham.id.au/paste/explain.pl?regex=",
        utils::URLencode(pattern))

    ## replace invalid characters
    chars <- c(";", "+", "&")
    reps <- c("%3B", "%2B", "%26")

    for (i in seq_along(reps)){
        URL2 <- gsub(chars[i], reps[i], URL2, fixed=TRUE)
    }

    lns <- try(suppressWarnings(readLines(URL2)), TRUE)

    if (length(lns) == 1 && grepl("Error in file", lns, TRUE)) {
        stop("Cound not parse `pattern`. Check your Internet connection.")
    }

    lns <- gsub("&quot;", "\"", lns[grep("NODE", lns):(length(lns) - 2)], fixed=TRUE)
    lns <- gsub("&gt;", ">", gsub("&lt;", "<", lns, fixed=TRUE), fixed=TRUE)
    lns <- gsub("\\", "\\\\", lns, fixed=TRUE)
    lets <- c("n", "r", "t", "f", "a")
    for (i in seq_len(length(lets))){
        lns <- gsub(paste0("\\\\", lets[i]), paste0("\\", lets[i]), lns, fixed=TRUE)
    }
    lns[length(lns)] <- gsub("</pre>$", "", lns[length(lns)])

    lns
}
