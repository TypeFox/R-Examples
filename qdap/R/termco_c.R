#' Combine Columns from a termco Object
#' 
#' Combines the columns of a termco object.  Generally intended for internal 
#' use but documented for completeness.
#' 
#' @param termco.object An object generated by either \code{\link[qdap]{termco}}, 
#' \code{\link[qdap]{termco_d}} or \code{\link[qdap]{termco_c}}.
#' @param combined.columns The names/indexes of the columns to be combined.
#' @param new.name A character vector of length one to name the new combined 
#' column.
#' @param short.term logical.  If \code{TRUE} column names are trimmed versions 
#' of the match list, otherwise the terms are wrapped with 'term(phrase)'
#' @param zero.replace Value to replace zeros with.
#' @param elim.old logical.  If \code{TRUE} eliminates the columns that are 
#' combined together by the named match.list.
#' @param percent logical.  If \code{TRUE} output given as percent.  If 
#' \code{FALSE} the 
#' output is proportion.
#' @param digits Integer; number of decimal places to round when printing.   
#' @return Returns a return a list, of class \code{"termco"}, of data frames and 
#' information regarding word counts:
#' \item{raw}{raw word counts by grouping variable} 
#' \item{prop}{proportional word counts by grouping variable; proportional to 
#' each individual's word use} 
#' \item{rnp}{a character combination data frame of raw and proportional}     
#' \item{zero_replace}{value to replace zeros with; mostly internal use}   
#' \item{percent}{The value of percent used for plotting purposes.} 
#' \item{digits}{integer value od number of digits to display; mostly internal 
#' use}  
#' @seealso \code{\link[qdap]{termco}}
#' @export
termco_c <-
function(termco.object, combined.columns, new.name, short.term = TRUE,
    zero.replace = NULL, elim.old = TRUE, percent = NULL, digits = 2) { 
    if (!class(termco.object) %in% c("termco")) {
        stop("termco.object must be a termco object")
    }
    subdf <- function(df, ii) {
        do.call("data.frame", c(as.list(df)[ii, drop=FALSE], check.names=FALSE))
    }
    if (is.null(percent))  {
        percent <- termco.object[["percent"]]
    }
    if (is.null(zero.replace)) {
        zero.replace <- termco.object$zero_replace
    }
    if (is.null(digits)) {
        digits <- termco.object$digits
    }
    x <- termco.object$raw
    nms <- gsub("term(", "", colnames(x), fixed = TRUE)
    lens <- sapply(nms[-c(1:2)], nchar)
    nms[-c(1:2)] <- substring(nms[-c(1:2)], 1, lens - 1)
    colnames(x) <- nms 
    x2 <- qcombine(x, combined.columns = combined.columns, elim.old = elim.old)
    y2 <- x2[, -c(1:2), drop = FALSE]/x[, 2]

    ## Added 1-21-14 to deal with missing values
    y2[] <- lapply(y2,  nan2zero)
    if (percent){
        y2 <- y2*100
    }
    y2 <- data.frame(x2[, 1:2], y2, check.names = FALSE)
    rnp <- raw_pro_comb(x2[, -c(1:2), drop = FALSE], 
        y2[, -c(1:2), drop = FALSE], digits = digits, 
        percent = percent, zero.replace = zero.replace, override = TRUE)  
    rnp <- data.frame(x2[, 1:2], rnp, check.names = FALSE) 
    o <- list(raw = x2, prop = y2, rnp = rnp, zero.replace = zero.replace,
        percent = percent, digits = termco.object$digits)
    if (!short.term) {
        nms <- colnames(o[["raw"]])
        nms[-c(1:2)] <- paste0("term(", nms[-c(1:2)], ")")
        o[1:3] <- lapply(o[1:3], function(x){
            colnames(x) <- nms
            x
        })
    }
    class(o) <- "termco"
    return(o)
}


