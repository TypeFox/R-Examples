#' Compute the lengths and values of runs of equal values in a vector
#'
#' Similar \code{\link{rle}} in the base package but for \code{\link{ff}} 
#' vectors.
#'
#' @param x an \code{\link{ff}} vector
#' @param ... further arguments are passed on the \code{\link{chunk}}
#'
#' @return An object of class \code{rle} which is a list with components
#' \describe{
#'   \item{lengths}{an integer vector containing the length of each run.}
#'   \item{values}{a vector of the same length as `lenghts' with the 
#'      corresponding values.}
#' }
#' @note The resulting rle object is a memory object and must fit into memory.
#'
#' @seealso \code{\link{rle}} for an implementation that runs on ordinary vectors.
#' 
#' @export
rle_ff <- function(x, ...) {

    #chunks <- chunk(x, ...)
    chunks <- chunk(x)
    rles   <- vector(length(chunks), mode="list")
    for (i in seq_along(chunks)) {
        rles[[i]] <- do.call(cbind, rle(x[chunks[[i]]]))
        # check for overlap with previous chunk
        if (i > 1 && rles[[i]][1,2] == rles[[i-1]][nrow(rles[[i-1]]),2]) {
            rles[[i-1]][nrow(rles[[i-1]]),1] <- rles[[i-1]][nrow(rles[[i-1]]),1] +
                rles[[i]][1,1]
            rles[[i]] <- rles[[i]][-1,]
        }
    }

    rles    <- do.call(rbind, rles)
    values  <- rles[,2]
    lengths <- rles[,1]
    structure( list(lengths=lengths, values=values) , class="rle")
}


# TODO export this stuff
rle <- function(x, ...){
  UseMethod('rle')
}
rle.default <- base::rle
rle.ff <- rle_ff
