# Hold-out method
# @param df A data-frame
hold_out <- function(df) {
    i <- 0
    imax <- length(df[,1])
    function() {
        i <<- i + 1
        if(i > imax) return(NA)
        one <- df[i,]
        rest <- df[-i,]
        return(list(one=one,rest=rest))
    }
}
