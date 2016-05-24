"year.strict" <- function(...)
    stop("you must expand 2-digit year abbreviations")

"year.expand" <-
function(y, cut.off = 30, century = c(1900, 2000), ...)
{
    ## cut.off specifies year for rounding up/down
    if(!is.numeric(y))
        stop("must be a numeric year specification")
    i <- (!is.na(y) & (y >= 0) & (y <= 99))
    if(any(i))
        y[i] <- ifelse(y[i] < cut.off,
                       y[i] + century[2],
                       y[i] + century[1])
    y
}
