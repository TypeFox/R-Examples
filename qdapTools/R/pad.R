#' Pad Strings
#' 
#' A convenience wrapper for \code{\link[base]{sprintf}} that pads the front end 
#' of strings with spaces or 0s. Useful for creating multiple uniform directories 
#' that will maintain correct order.
#' 
#' @param x A character, factor, numeric vector.
#' @param padding Number of characters to pad.  Default makes all elements of a 
#' string the number of characters of the element with the maximum characters.
#' @param sort logical.  If \code{TRUE} the outcome is sorted.
#' @param type A character string of \code{"detect"}, \code{"numeric"}, 
#' \code{"character"}, \code{"d"} or \code{"s"}.  If numeric zeros are padded.  
#' If character spaces are padded.  The \code{detect} attempts to determine if x 
#' is numeric (d) or not (s).
#' @return Returns a character vector every element padded with 0/spaces.
#' @note \code{pad} is a wrapper for the \code{\link[base]{sprintf}} function.
#' \code{pad} may behave differently on various platforms in accordance with the
#' documentation for \code{\link[base]{sprintf}}: "actual implementation will 
#' follow the C99 standard and fine details (especially the behaviour under user 
#' error) may depend on the platform."  See \code{\link[base]{sprintf}} for more 
#' information.
#' @export
#' @seealso \code{\link[base]{sprintf}}
#' @examples
#' pad(sample(1:10, 10))
#' pad(sample(1:10, 10), sort=FALSE)
#' pad(as.character(sample(1:10, 10)))
#' pad(as.character(sample(1:10, 10)), sort=FALSE)
#' pad(as.character(sample(1:10, 10)), 4)
#' pad(month.name)
pad <- function (x, padding = max(nchar(as.character(x))), sort = TRUE, 
    type = "detect") {
    poss <- c("detect", "numeric", "character", "d", "s")
    if (!type %in% poss) 
        stop("type must be: \"detect\", \"numeric\"\\\"d\" or \"character\"\\\"s\"")
    Rel <- c(NA, "d", "s", "d", "s")
    type <- Rel[poss %in% type]
    if (is.na(type)) {
        type <- ifelse(is.numeric(x), "d", "s")
    }
    x <- sprintf_ish(x, padding, type)
    if (sort) {
        x <- sort(x)
    }
    x
}

sprintf_ish <- function(x, padding, type){
    OS <- Sys.info()[['sysname']]
 
    if (OS %in% c("Windows", "Linux")) {
        sprintf(paste0("%0", padding, type), x)
    } else {
        type <- ifelse(type == "s", " ", "0")
        pads <- sapply(padding - nchar(x), function(i)  {
            if (i == 0) return("")
            paste(rep(type, i), collapse = "")
        })
        paste0(pads, x)
      
    }
}




