# Verified 1.3.18
# Version 3.0
summary.Catalog <-
function(object, ...) {
        b = object$catalog
        if ( class(b[[1]]) != "list" ) { b = list(); b[[1]] = object$catalog }
        z = plyr::ldply(b, function(x) { c(
                x$Name, 
                x$State, 
                round(x$Longitude, 2), 
                round(x$Latitude, 2), 
                x$Altitude, 
                paste0(x$Install, collapse="/"), 
                paste0(x$Start, collapse="/"), 
                paste0(paste0(range(x$Avble.yrs), collapse=" -> "), " (", length(x$Avble.yrs), ")"), 
                paste0(x$Measure.code, ":", x$Measure.unit)) })
        std <- c("Name","Altitude","Latitude","Longitude","Measure.unit","Measure.code","Install","Start","State","Avble.yrs")  
        colnames(z) <- c("Name", "State", "Longitude", "Latitude", "Altitude", "Inst", "Start","Data range (yrs)", "Measurements")
        extra = !( names(b[[1]]) %in% std )
        z.ex = character(0)
        if ( any(extra) ) {
                extra = names(b[[1]])[extra]
                z.ex = apply(array(extra), 1, function(y, b) { plyr::ldply(b, function(x, y) { get(y, x) }, y) }, b)
                z.ex = t(plyr::ldply(z.ex, function(x) { t(x) } ))
                colnames(z.ex) <- extra
                z = cbind(z, z.ex)
        }
        (z)
}
