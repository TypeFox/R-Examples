# Version 1.3.21
# CC3
Catalog <- function(catalog, data, ...){
        if ( any(unlist(lapply(catalog, .canbeCatalog)) == FALSE) ) { stop('catalog cannot be coerced to class Catalog.') }
        if ( length(catalog) != length(data) ) { stop('catalog and data length differ.') }
        col <- list(catalog = catalog, data = data, ...)
        class(col) <- 'Catalog'
        return(col)
}

.canbeCatalog <- function(catalog) {
        NAMES <- names(catalog)
        is = FALSE
        if (    ( 'Name' %in% NAMES ) &
                ( 'Altitude' %in% NAMES ) &
                ( 'Latitude' %in% NAMES ) &
                ( 'Longitude' %in% NAMES ) &
                ( 'Measure.code' %in% NAMES ) &
                ( 'Measure.unit' %in% NAMES ) &
                ( 'Install' %in% NAMES ) &
                ( 'Start' %in% NAMES ) &
                ( 'State' %in% NAMES ) &
                ( 'Avble.yrs' %in% NAMES )
        ) { is = TRUE }
        return(is)
}
