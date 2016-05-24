summary.hwtANYN <-
function (object, ...) 
{
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Filter was: Haar\n")
    cat("Transform type: ", object$type, "\n")
    cat("Object was reindex to match wd: ", object$reindex, "\n")
}
