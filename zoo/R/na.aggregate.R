na.aggregate <- function(object, ...) UseMethod("na.aggregate")

## fills NA values with some aggregated function of the data.
## generalises imputing by the overall mean, by calendar month, etc.
na.aggregate.default <- function(object, by = 1, ..., FUN = mean, na.rm = FALSE, maxgap = Inf)
{
    if (is.function(by)) by <- by(time(object), ...)
    ## applied to each aggregated group in each series:
    f <- function(x)
        replace(x, is.na(x), FUN(x[!is.na(x)]))
    na.aggregate.0 <- function(y) {
        yf <- ave(y, by, FUN = f)
        .fill_short_gaps(y, yf, maxgap = maxgap)
    }
    object[] <- if (length(dim(object)) == 0) na.aggregate.0(object)
                else apply(object, 2, na.aggregate.0)
    if (na.rm) na.trim(object, is.na = "all") else object
}
