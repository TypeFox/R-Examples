residuals.bayesx <- function(object, model = NULL, term = NULL, ...)
{
  object <- get.model(object, model)
  rval <- list()
  k <- length(object)
  mn <- rep("model", length.out = k)
  for(i in 1L:k) {
    if(!is.null(term)) {
      if(length(term) > 1L) {
        rval[[i]] <- list()
        tn <- names(object[[i]]$effects)
        if(is.null(tn))
          tn <- paste(term)
        for(j in 1:length(term)) {
          tmp <- attr(object[[i]]$effects[[term[j]]], "partial.resids")
          eval(parse(text = paste("rval[[i]]$'", tn[j], "' <- tmp", sep = "")))
        }
      } else rval[[i]] <- attr(object[[i]]$effects[[term]], "partial.resids")
    } else {
      rval[[i]] <- bayesx.reorder(object[[i]], object[[i]]$residuals, TRUE)
    }
    if(!is.null(object[[i]]$bayesx.setup$model.name))
      mn[i] <- object[[i]]$bayesx.setup$model.name
  }
  mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
  names(rval) <- mn
  if(length(rval) < 2L)
    rval <- rval[[1L]]
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA
  if(any(is.na(rval)))
    warning("residuals are missing in object!")

  return(rval)
}

