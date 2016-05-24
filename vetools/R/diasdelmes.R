# Verified 1.3.18
diasdelmes <-
function(y, meses) {
        if ( any( meses > 23 ) ) { stop("meses can not be greater than 23 months.") }
        x <- array(data = meses)
        z <- apply(x, 1, function(m, y){
                if (m >= 12) { y = y + 1; m = m - 12; if(m == 0) {m = 1} }
                as.POSIXlt(seq(as.Date(paste(y, m+1, "01",sep="-")), by="month", length.out=1)-1)$mday
        }, y = y)
        return(sum(z))
}
