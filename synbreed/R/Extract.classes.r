"[.relationshipMatrix" <- function(x,...) {
   y <- NextMethod("[")
   class(y) <- oldClass(x)
   y
}

"[.GenMap" <- function(x,...) {
   y <- NextMethod("[")
   class(y) <- oldClass(x)
   y
}
