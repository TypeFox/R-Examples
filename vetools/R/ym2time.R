# Verified 1.3.18
ym2time <-
function(e) {
        if (length(e) != 2) { 
                stop("Wrong size: must be vector of length two, c(year, month).")
        }
        return( tis::jul(e[1] + ((e[2]-1)/12)) )  
}
