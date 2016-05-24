impute <- function(x, what=c("median", "mean")){

    what <- match.arg(what)
    
    if(what == "median"){
        retval <-
            apply(x, 2,
                  function(z) {z[is.na(z)] <- median(z, na.rm=TRUE); z})
    }
    else if(what == "mean"){
        retval <-
            apply(x, 2,
                  function(z) {z[is.na(z)] <- mean(z, na.rm=TRUE); z})
    }
    retval
}
