compact <-
function(x) { 
    null <- vapply(x, is.null, logical(1)) 
    x[!null] 
}
