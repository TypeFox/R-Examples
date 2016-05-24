# Verified 1.3.18
time2ym <-
function(d) {
        if ( !("Date" %in% class(d)) ) { stop("Not a Date class object") }
        return(c(as.numeric(format(d, format=c("%Y"))), as.numeric(format(d, format=c("%m")))))
}
