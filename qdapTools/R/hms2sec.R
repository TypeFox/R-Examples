#' Convert h:m:s To/From Seconds
#' 
#' \code{hms2sec} - Converts a vector of h:m:s to seconds.
#' 
#' @param x A vector of times in h:m:s (for \code{hms2sec}) or seconds (for 
#' \code{sec2hms}) .
#' @return \code{hms2sec} - Returns a vector of times in seconds.  
#' @keywords time, conversion
#' @seealso \code{\link[chron]{times}}
#' @rdname timeconv
#' @export
#' @examples 
#' hms2sec(c("02:00:03", "04:03:01"))
#' hms2sec(sec2hms(c(222, 1234, 55)))
#' sec2hms(c(256, 3456, 56565))
hms2sec <- 
function(x) {
    hms <- as.character(x)
    op <- FALSE
    if (length(hms) == 1) {
        hms <- c(hms, "00:00:00")
        op <- TRUE  
    }
    DF <- sapply(data.frame(do.call(rbind, strsplit(hms, ":"))), function(x){
        as.numeric(as.character(x))
    })
    out <- DF[, 1] * 3600 + DF[, 2] * 60 + DF[, 3]
    if (op) {
        out <- out[1]
    }
    out
}



#' Convert Seconds to h:m:s
#' 
#' \code{sec2hms} - Converts a vector of seconds to h:m:s.
#' 
#' @return \code{sec2hms} - Returns a vector of times in h:m:s format.  
#' @export
#' @rdname timeconv
#' @importFrom chron times
sec2hms <-
function(x) {
    l1 <- FALSE
    if (length(x) == 1) {
        x <- c(x, 0)
        l1 <- TRUE
    } 
    h <- floor(x/3600)
    m <- floor((x-h*3600)/60)
    s <- x-(m*60 + h*3600)
    pad <- function(x) sprintf("%02d", as.numeric(x))
    out <- times(paste2(data.frame(apply(data.frame(h=h, m=m, s=s), 
       2, pad)), sep=":"))
    if (l1) {
        out <- out[1]
    }
    out
}


