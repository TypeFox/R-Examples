is.na.matrix <-
function (mat, rows = TRUE) 
{
    if (!is.matrix(mat)) 
        stop("'mat' must be a matrix or data frame.")
    if (rows) 
        return(apply(mat, 1, function(x) any(is.na(x))))
    else return(apply(mat, 2, function(x) any(is.na(x))))
}
