is.Catalog <- function(x, ignore.class = FALSE) {
        if ( any(unlist(lapply(x$catalog, .canbeCatalog)) == FALSE) == TRUE ) { return(FALSE) }
        if ( length(x$catalog) != length(x$data) ) { return(FALSE) }
        if ( ( ! ignore.class ) & ( class(x) != 'Catalog' ) ) { return(FALSE) }
        return(TRUE)
}