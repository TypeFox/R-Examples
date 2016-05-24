as.Catalog <- function(x) {
        if ( is.Catalog(x, ignore.class = TRUE) == FALSE ) { stop('Cannot coherce x to class Catalog.') }
        class(x) <- 'Catalog'
        return(x)
}