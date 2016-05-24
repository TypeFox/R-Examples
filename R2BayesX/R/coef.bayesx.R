coef.bayesx <- function(object, model = NULL, ...)
{
  object <- get.model(object, model)
  if(length(object) > 1L) {
    rval <- vector("list", length = length(object))
    for(i in 1L:length(object))
      rval[[i]] <- object[[i]]$fixed.effects
    names(rval) <- names(object)
  } else rval <- object[[1L]]$fixed.effects
  cattr <- attributes(rval)
  rvalo <- rval
  rval <- as.data.frame(rval)
  eattrn <- names(cattr)
  for(i in 1L:length(eattrn))
    if(eattrn[i] != "dim" && eattrn[i] != "dimnames" && eattrn[i] != "names"
      && eattrn[i] != "row.names" && eattrn[i] != "class")
      attr(rval, eattrn[i]) <- attr(rvalo, eattrn[i])

  return(rval)
}

