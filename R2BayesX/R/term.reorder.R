term.reorder <- function(x, info) 
{
  if(!is.null(x) && length(x) > 0L) {
    rval <- list()
    info <- readLines(info)
    ni <- length(info) - 1L
    nx <- length(x)
    nt <- rep(NA, ni)
    taken <- NULL
    for(k in 1L:ni) {
      pt <- try(parse(text = info[k]), silent = TRUE)
      if(inherits(pt, "try-error")) next
      term <- eval(pt)
      if(is.null(term$term)) next
      if(!is.null(term$by) && term$by != "NA")
        term$term <- paste(term$term, ":", term$by, sep = "")
      if(is.null(term$isFactorBy))
        term$isFactorBy <- FALSE
      if(is.rt(term$term))
        term$term <- rrd(term$term)
      nt[term$pos] <- term$term
      if(!term$isFactor && !term$isFactorBy) {
        for(j in 1L:nx) {
          sl2 <- splitme(attr(x[[j]], "specs")$label)
          if(length(sl2) > 1L)
            label2 <- paste("te", resplit(sl2[2L:length(sl2)]), sep = "") 
          else
            label2 <- "not.in.labels"
          if((term$term == attr(x[[j]], "specs")$label && term$class %in% class(x[[j]])) || term$term == label2) {
            if(!is.null(term$map)) {
              if(length(grep(term$map, ls(envir = globalenv()))))
                attr(x[[j]], "map.name") <- term$map
            }
            rval[[term$pos]] <- x[[j]]
            taken <- c(taken, j)
          }
        }
      } 
      if(term$isFactor && !term$isFactorBy) {
        fc <- list()
        idf <- 1L
        if(is.null(term$realname))
          tnames <- rrmfs(term$names)
        else
          tnames <- term$realname
        for(j in 1L:nx) {
          if(attr(x[[j]], "specs")$label %in% tnames) {
            if(!is.null(term$map)) {
              if(length(grep(term$map, ls(envir = globalenv()))))
                attr(x[[j]], "map.name") <- term$map
            }
            fc[[idf]] <- x[[j]]
            idf <- idf + 1L
            taken <- c(taken, j)
          }
        }
        attr(fc, "specs") <- list(dim = 1L, term = term$term, 
          label = term$term, is.factor = TRUE)
        class(fc) <- "linear.bayesx"
        rval[[term$pos]] <- fc
      }
      if(!term$isFactor && term$isFactorBy) {
        by <- list()
        xl <- names(x)
        for(k in 1L:length(term$isFactorByNames))
          if(!is.na(j <- match(term$isFactorByNames[k], xl))) {
            if(!is.null(term$map)) {
              if(grep(term$map, ls(envir = globalenv())))
                attr(x[[j]], "map.name") <- term$map
            }
            by[[xl[j]]] <- x[[j]]
            taken <- c(taken, j)
          }
        class(by) <- "varying.bayesx"
        rval[[term$pos]] <- by
      }
    }
    if(length(rval) > 0L)
      try(names(rval) <- nt, silent = TRUE)
    if(length(taken) < nx) {
      nterms <- names(x)[!c(1L:nx) %in% taken]
      nrval <- list()
      ntn <- NULL
      for(k in 1L:length(nterms)) {
        nrval[[k]] <- x[[nterms[k]]]
        ntn <- c(ntn, nterms[k])
      }
      try(names(nrval) <- ntn, silent = TRUE)
      rval <- c(rval, nrval)
    }
    for(k in 1:length(rval)) {
      if(!is.null(rval[[k]]) && !is.list(rval[[k]]) && nrow(rval[[k]]) < 2L) 
        attr(rval[[k]], "specs")$is.factor <- TRUE
    }
    namrval <- names(rval)
    rval <- x2df(rval, rn = TRUE)
    if(!inherits(rval, "list")) {
      rval <- list(rval)
      names(rval) <- namrval
    }
    for(k in 1:length(rval)) {
      if(length(attr(rval[[k]], "specs")$label))
        names(rval)[k] <- attr(rval[[k]], "specs")$label
    }
    return(rval)
  } else {
    xnames <- names(x)
    x <- x2df(x, rn = TRUE)
    if(!inherits(x, "list")) {
      x <- list(x)
      names(x) <- xnames
    }
    return(x)
  }
}

