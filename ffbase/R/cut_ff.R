#' cut divides the range of x into intervals and codes the values in x according to which interval they fall. The leftmost interval corresponds to level one, the next leftmost to level two and so on. 
#'
#' The \code{cut} method for ff with the behaviour of \code{link{cut}}
#' @title Convert Numeric ff vector to factor ff
#' @export
#' @export cut.ff
#' @seealso cut
#' @method cut ff
#' 
#' @param x a (numeric) ff object that will be cut into pieces
#' @param breaks specifies the breaks for cutting this
#' @param ... other parameters that can be given to \code{\link{cut.default}}
#' 
#' @return ff a new \code{\link{ff}} object with the newly created factor
cut.ff <- function(x, breaks, ...){
   f <- NULL
   
   #### borrowed code from cut.default
   if (length(breaks) == 1L) {
     if (is.na(breaks) || breaks < 2L) 
       stop("invalid number of intervals")
     nb <- as.integer(breaks + 1)
     dx <- diff(rx <- range(x, na.rm = TRUE))
     if (dx == 0) {
       dx <- abs(rx[1L])
       breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
                         length.out = nb)
     }
     else {
       breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
       breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
                                dx/1000)
     }
   }
   else nb <- length(breaks <- sort.int(as.double(breaks)))
   if (anyDuplicated(breaks)) 
     stop("'breaks' are not unique")
   ####
   
   args <- list(...)
   args$by <- NULL
   args$breaks <- breaks
   for (i in chunk(x, ...)){
     Log$chunk(i)
     args$x <- x[i]
     res <- do.call(cut, args)
     f <- ffappend( f
                  , res
                  , adjustvmode=FALSE
                  )
   }
   f
}
