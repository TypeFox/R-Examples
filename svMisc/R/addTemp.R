addTemp <- function (x, item, value, use.names = TRUE, replace = TRUE)
{
    x <- as.character(x)[1]
    item <- as.character(item)[1]
    if (existsTemp(x)) dat <- getTemp(x) else dat <- list()
    ## The object must be a list!
    if (!inherits(dat, "list")) stop(x, " must be a list!")
    ## Does 'item' already exist?
	if (item %in% names(dat))
	    value <- addItems(dat[[item]], value,
			use.names = use.names, replace = replace)
    dat[[item]] <- value
	assignTemp(x, dat)
}
