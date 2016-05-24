changeTemp <- function (x, item, value, replace.existing = TRUE)
{
    x <- as.character(x)[1]
    item <- as.character(item)[1]
    if (existsTemp(x)) dat <- getTemp(x) else dat <- list()
    ## The object must be a list!
    if (!inherits(dat, "list")) stop(x, " must be a list!")
    ## Does 'item' already exist?
    if (replace.existing || !item %in% names(dat)){
        dat[[item]] <- value
		assignTemp(x, dat)
    }
}
