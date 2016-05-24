"fun.zero.omit" <-
function(object)
{
d <- dim(object)
if(length(d) != 2) {
object <- as.matrix(object)
d <- dim(object)
}
if(length(nas <- fun.which.zero(object))) {
nas <- unique((nas - 1) %% d[1] + 1)
object <- object[ - nas,  , drop = FALSE]
}
attr(object, "na.message") <- paste("Dropped", length(nas), "cases due to zero values")
object
}

