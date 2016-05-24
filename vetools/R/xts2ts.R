# Verified 1.3.18
xts2ts <-
function(b.xts) {
        b.xts.na = b.xts[!is.na(b.xts)]
        datos = c()
        datos[1] = b.xts.na[1]
        k = 2
        last = 1
        extra = c()
        for ( i in 2:length(b.xts.na) ) {
                d = diffmonths(b.xts.na[i], b.xts.na[i-1])              
                if ( d > 1 ) {
                        datos[k + d - 1] <- b.xts.na[i]
                        last = k + d - 1
                        k <- k + d
                } else if( d == 1 ) {
                        datos[k] <- b.xts.na[i]
                        last = k
                        k = k + 1
                } else {
                        if (datos[last] != b.xts.na[i]) {
                                extra = c(extra, b.xts.na[i])
                                warning("Two or measurements in same date do not match, averaging measurements.")
                                datos[last] = 0.5 * ( datos[last] + b.xts.na[i] )
                        }
                }
        }
        b.ts = ts(datos, start=time2ym(start(b.xts.na)), frequency=12)
        
        if ( any(start(b.ts) != c(as.numeric(format(start(b.xts), "%Y")), as.numeric(format(start(b.xts), "%m")))) ) { stop("PANIC: Start date conversion did not succeed.") }
        if ( any(end(b.ts) != c(as.numeric(format(end(b.xts), "%Y")), as.numeric(format(end(b.xts), "%m")))) ) {stop("PANIC: End date conversion did not succeed.") }
        return(b.ts)
}
