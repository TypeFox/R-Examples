# -----------------------------------------------------
# This function displays the size of the objects
# Allows to detect the memory consuming objects
# -----------------------------------------------------
all.object.sizes <- function() {
  res <- NULL
  nomres <- NULL
  for (x in ls(envir=parent.frame())) {
     z<- mget(x, envir=parent.frame(), ifnotfound = list(NULL))
    if (!is.null(zu<-unlist(z, use.names = FALSE))) {
#     res <- c(res, object.size(mget(x, envir=parent.frame(), ifnotfound = list(1))))
     res <- c(res, object.size(zu))
     nomres <- c(nomres, x)
     } else cat(x, "not found trouve\n")
   } # end for
  names(res) <-nomres
  return(sort(res))
}
