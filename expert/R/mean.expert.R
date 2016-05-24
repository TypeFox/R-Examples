### ===== expert =====
###
### Mean of the aggregated expert distribution
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>
###          Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

mean.expert <- function(x, ...)
{
    ## Get bounds
    breaks <- x$breaks

    ## Compute group midpoints
    midpoints <- breaks[-length(breaks)] + diff(breaks)/2

    ## Drop the boundaries column and convert to matrix for use in
    ## crossprod()
    x <- as.matrix(x$probs)

    ## Compute mean per column
    drop(crossprod(x, midpoints))
}
