fitted.bayesx <- function(object, model = NULL, term = NULL, ...)
{  
  object <- get.model(object, model)
  k <- length(object)
  mn <- rep("model", length.out = k)
  if(is.null(term)) {
    if(length(object) > 1L) {
      rval <- vector("list", length = k)
      for(i in 1L:k) {
        rval[[i]] <- bayesx.reorder(object[[i]], object[[i]]$fitted.values, TRUE)
        if(!is.null(object[[i]]$bayesx.setup$model.name))
          mn[i] <- object[[i]]$bayesx.setup$model.name
      }
      mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
      names(rval) <- mn
    } else {
      rval <- bayesx.reorder(object[[1L]], object[[1L]]$fitted.values, TRUE)
    }
  } else {
    if(length(object) > 1L) {
      rval <- list()
      for(i in 1L:k) {
        rval[[i]] <- object[[i]]$effects[term]
        if(!is.null(object[[i]]$bayesx.setup$model.name))
          mn[i] <- object[[i]]$bayesx.setup$model.name
      }
      mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
      names(rval) <- mn
    } else {
      if(length(term) > 1L) {
        rval <- object[[1L]]$effects[term]
      } else {
        rval <- object[[1L]]$effects[[term]]
      }
    }
  }
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA
  else {
    if(is.list(rval) && any(grepl("bayesx", class(rval[[1L]]))) || 
      (is.list(rval) && is.list(rval[[1L]])))
      class(rval) <- "fit.bayesx"
  }
  if(any(is.na(rval)))
    warning("fitted values are missing in object!")
  rval <- x2df(rval)

  return(rval)
}


"[.fit.bayesx" <- function(x, term)
{
  if(is.list(x)) {
    if(is.character(term))
      if(any(is.na(term <- pmatch(term, names(x)))))
        stop("element not existing!")
    return(x[[term]])
  } else return(x)
}


x2df <- function(x, rn = FALSE)
{
  if(!is.data.frame(x) && !is.null(x) && length(x)) {
    if(is.list(x)) {
      for(i in 1L:length(x)) {
        if(is.list(x[[i]])) {
          x[[i]] <- x2df(x[[i]])
        } else {
          if(!is.data.frame(x[[i]]) && !is.null(x[[i]]) & !is.null(dim(x[[i]]))) {
            xattr <- attributes(x[[i]])
            nxa <- names(xattr)
            cx <- class(x[[i]])
            if(any(grepl("bayesx", cx, fixed = TRUE)))
              cx <- grep("bayesx", cx, fixed = TRUE, value = TRUE)
            if(rn)
              x[[i]] <- chacol(x[[i]])
            x[[i]] <- as.data.frame(x[[i]])
            class(x[[i]]) <- c(cx, class(x[[i]]))
            for(k in 1L:length(nxa))
              if(all(nxa[k] != c("dim", "dimnames", "class", "names", "row.names")))
                attr(x[[i]], nxa[k]) <- xattr[[k]]
          }
        } 
      }
    if(length(x) < 2L) {
      if(is.null(attr(x, "specs")$is.factor))
        x <- x[[1L]]
    }
    } else {
      if(!is.null(dim(x))) {
        xattr <- attributes(x)
        nxa <- names(xattr)
        cx <- class(x)
        if(any(grepl("bayesx", cx, fixed = TRUE)))
          cx <- grep("bayesx", cx, fixed = TRUE, value = TRUE)
        if(rn)
          x <- chacol(x)
        x <- as.data.frame(x)
        class(x) <- c(cx, class(x))
        for(k in 1L:length(nxa)) 
          if(all(nxa[k] != c("dim", "dimnames", "class", "names", "row.names")))
            attr(x, nxa[k]) <- xattr[[k]]
      }
    }
  }

  return(x)
}


bayesx.reorder <- function(object = NULL, x, unique = FALSE) {
  if(!is.null(object)) {
    i <- if(is.list(object)) object$bayesx.setup$order else object
    if(!is.null(i)) {
      if(!is.null(dim(x))) {
        j <- 1:nrow(x)
        x <- x[j[order(i)], ]
        rownames(x) <- 1:nrow(x)
      } else {
        j <- 1:length(x)
        x <- x[j[order(i)]]
      }
    }
    if(is.list(object)) {
      if(!is.null(object$bayesx.setup$YLevels)) {
        if(object$bayesx.setup$Yn %in% colnames(x)) {
          x[[object$bayesx.setup$Yn]] <- factor(x[[object$bayesx.setup$Yn]],
            levels = as.integer(object$bayesx.setup$nYLevels),
            labels = object$bayesx.setup$YLevels)
        }
      }
    }
  }
  if(unique & !is.null(dim(x))) {
    if(ncol(x) == 2) {
      ui <- apply(x, 1, function(x) any(duplicated(x)))
      if(all(ui))
        x <- as.numeric(unlist(x[, 1L]))
    }
  }

  return(x)
}

