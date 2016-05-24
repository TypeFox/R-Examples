#' Tabulate Frequency Counts for Multiple Vectors 
#' 
#' Similar to \code{\link[base]{tabulate}} that works on multiple vectors.
#' 
#' @param vects A \code{\link[base]{vector}}, \code{\link[base]{list}}, or 
#' \code{\link[base]{data.frame}} of named/unnamed vectors.
#' @keywords tabulate frequency
#' @export
#' @seealso \code{\link[base]{tabulate}}, \code{\link[qdapTools]{counts2list}}
#' @return Returns a \code{\link[base]{data.frame}} with columns equal to 
#' number of unique elements and the number of rows equal to the the original 
#' length of the \code{\link[base]{vector}}, \code{\link[base]{list}}, or 
#' \code{\link[base]{data.frame}} (length equals ncols in 
#' \code{\link[base]{data.frame}}).  If list of vectors is named 
#' these will be the rownames of the dataframe.
#' @author Joran Elias and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @references \url{http://stackoverflow.com/a/9961324/1000343}
#' @examples 
#' mtabulate(list(w=letters[1:10], x=letters[1:5], z=letters))
#' mtabulate(list(mtcars$cyl[1:10]))
#' 
#' ## Dummy coding
#' mtabulate(mtcars$cyl[1:10])
#' mtabulate(CO2[, "Plant"])
#' 
#' dat <- data.frame(matrix(sample(c("A", "B"), 30, TRUE), ncol=3))
#' mtabulate(dat)
#' t(mtabulate(dat))
#' counts2list(mtabulate(dat))
mtabulate <- function(vects) { 
    lev <- sort(unique(unlist(vects)))
    dat <- do.call(rbind, lapply(vects, function(x, lev){ 
        tabulate(factor(x, levels = lev, ordered = TRUE),
        nbins = length(lev))}, lev = lev))
    colnames(dat) <- sort(lev) 
    data.frame(dat, check.names = FALSE)
}


## param drop.na logical.  If \code{TRUE} \code{NA} columns (elements that 
## contained just an \code{NA}) will be dropped.## @author akrun of StackOverflow and Tyler Rinker <tyler.rinker@@gmail.com>.
## references url{http://stackoverflow.com/a/32753233/1000343}
## mtabulate <- function(vects, drop.na = TRUE) {
## 
##     x <- y <- . <- NULL
##     vects <- as.list(vects)
##     if (is.null(names(vects))) names(vects) <- seq_along(vects)
##     dat <- data.table::data.table(
##         x = names(vects),
##         y = vects, 
##         stringsAsFactors = FALSE
##     )
##     dat$y <- relist(unlist(dat$y), skeleton=dat$y)
##     data.table::setDT(dat)
##     dat <- dat[, .(y = unlist(y)), by = x]
##     out <- suppressMessages(data.table::dcast(dat, x ~ y, fun=length, drop=FALSE, fill=0))
##     out2 <- as.data.frame(out[, -1, with=FALSE])
##     rownames(out2) <- out[[1]]
##     if (isTRUE(drop.na)) out2[, "NA"] <- NULL
##     out2[names(vects), ]
## }
## }

