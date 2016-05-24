# Verified 1.3.18
est.fill <-
function(collection, cut = c(1968,3), at.start=T) {
        serie = collection$data
        for ( i in 1 : length(serie) ) {
                if ( at.start == T ) {
                        f1 = start(serie[[i]])
                } else {
                        f1 = end(serie[[i]])
                }
                f1 = f1[1] + f1[2]/12
                f2 = cut[1] + cut[2]/12
                if ( (at.start == T) &&  ( f1 > f2 ) ) {
                        window(serie[[i]], start=cut, end=start(serie[[i]])) = NA
                        serie[[i]] = fill.small.missing(serie[[i]], max.len=1e6)
                }
                
                if ( (at.start == F) &&  ( f1 < f2 ) ) {
                        window(serie[[i]], end=cut, start=end(serie[[i]]), extend=T) = NA
                        serie[[i]] = fill.small.missing(serie[[i]], max.len=1e6)
                }      
        }
        col <- list(data = serie, catalog = collection$catalog)
        class(col) <- "Catalog"
        return(col)
}
