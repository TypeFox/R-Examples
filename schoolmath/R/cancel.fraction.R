cancel.fraction <-
function (numerator, denominator) 
{
    test <- is.whole(numerator)
    if (test == FALSE) {
        msg <- cat("Please enter a whole numerator\r")
        return(msg)
    }
    test <- is.whole(denominator)
    if (test == FALSE) {
        msg <- cat("Please enter a whole denominator\r")
        return(msg)
    }
## The following line was modified.
    test <- is.prim(denominator)  ## Initially, is.prim(numerator), actually, I am not sure this test is important!
## End of modifications
    if (test == FALSE) {
        ggT <- gcd(numerator, denominator)
        numerator <- numerator/ggT
        denominator <- denominator/ggT
    }
    cat(numerator, "/", denominator)
}

