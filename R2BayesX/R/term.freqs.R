term.freqs <- function(object, model = NULL, term = NULL, ...)
{
  if(is.null(object) || length(object) < 2L)
    stop("nothing to do!")
  is.fit <- FALSE
  if(inherits(object, "data.frame")) {
    if(is.null(xn <- attr(object, "specs")$term))
      xn <- paste(xn, collapse = ".", sep = "")
    else xn <- names(object)[1L]
    object <- list(object)
    names(object) <- xn
    object <- list(list(effects = object))
    term <- 1L
    is.fit <- TRUE
  }
  if(inherits(object, "fit.bayesx") && !is.fit) {
    is.fit <- TRUE
    if(any(grepl("bayesx", class(object[[1L]]))))
      object <- list(object)
    k <- length(object)
    tmp <- list()
    for(i in 1L:k)
      tmp[[i]] <- list(effects = object[[i]])
    object <- tmp
    if(is.null(term))
      term <- 1L
  }
  if(!is.fit)
    object <- get.model(object, model)
  rval <- list()
  k <- length(object)
  mn <- rep("model", length.out = k)
  if(is.null(term))
    term <- names(object[[1]]$effects)
  for(i in 1L:k) {
    rval[[i]] <- list()
    for(j in 1:length(term)) {
      if(is.character(term[j]) & !is.null(neff <- names(object[[i]]$effects))) {
        term[j] <- gsub("[[:space:]]", "", term[j])
        term[j] <- neff[pmatch(term[j], neff)]
        if(is.na(term[j])) stop(paste("term", term[j], "does not exist, no frequencies available!"))
        tx <- term[j]
      } else {
        if(term[j] > length(object[[i]]$effects))
          stop("term does not exist, no frequencies available!")
        tx <- names(object[[i]]$effects)[as.integer(term[j])]
      }
      rval[[i]][[tx]] <- df <- attr(object[[i]]$effects[[tx]], "df")
    }
    if(!is.null(object[[i]]$bayesx.setup$model.name) && k > 1L)
      mn[i] <- object[[i]]$bayesx.setup$model.name
    if(length(term) < 2)
      rval[[i]] <- rval[[i]][[1L]]
  }
  if(k > 1L) {
    mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
    names(rval) <- if(is.null(names(object))) mn else names(object)
  } else rval <- rval[[1L]]
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA
  if(any(is.na(rval)))
    warning("frequency tables are missing in object!")
  rval <- delete.NULLs(rval)

  return(rval)
}

