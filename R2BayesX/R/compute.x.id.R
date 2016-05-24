compute.x.id <- function(x, id = NULL, c.select = NULL, range = NULL, symmetric = TRUE)
{ 
  if(is.null(id) && (is.vector(x) || is.array(x))) {
    if(!is.null(names(x))) {
      id <- names(x)
      x <- as.vector(x)
    }
  }
  if(is.factor(id))
    id <- as.character(id)
  if(is.array(x) && length(dim(x)) < 2L)
    x <- as.vector(x)
  if(is.null(dim(x)) && is.null(dim(id))) {
    if(length(x) != length(id))
      stop("arguments x and id are differing!")
  } else {
    x <- unclass(x)
    if(is.list(x)) 
      nx <- names(x)
    if(is.matrix(x)) {
      if(ncol(x) < 2 & !is.null(id)) {
        x <- data.frame("id" = id, "x" = as.numeric(x))
        nx <- names(x)
        c.select <- "x"
        id <- NULL
      } else {
        x <- as.list(as.data.frame(x))
        nx <- names(x)
        if(all(nx %in% paste("V", 1L:length(nx), sep = ""))) {
          nx[1L:2L] <- c("id", "x")
          c.select <- "x"
        }
      }
    }
    if(is.data.frame(x)) {
      x <- as.list(x)
      nx <- names(x)
    }
    if(is.null(id))
      id <- x[[1L]]
    else {
      if(is.character(id)) {
        if(is.na(id <- pmatch(id, nx)))
          stop("argument id is specified wrong!")
      } else {
        if(id > length(nx))
          stop("argument id is specified wrong!")
      }
      id <- x[[id]]
    }
    if(is.null(c.select)) {
      take <- c("mean", "Mean", "MEAN", "estimate", 
        "Estimate", "ESTIMATE", "mean", "pmode", "pmean_tot")
      did.take <- FALSE
      for(k in take) {
        if(!is.na(pmatch(k, nx)) & !did.take) {
          x <- x[[k]]
          did.take <- TRUE
        }
     }
     if(!did.take && length(x) > 1L)
       x <- x[[2L]]
    } else {
      if(is.character(c.select)) {
        k <- pmatch(c.select, nx)
      if(is.na(k))
        stop("argument c.select is specified wrong!")
      x <- x[[k]]
      } else {
        if(c.select > length(nx))
          stop("argument c.select is specified wrong!")
        x <- x[[c.select]]
      }
    }
  }
  xrange <- range(x, na.rm = TRUE)
  if(symmetric) {
    xrange <- c(-1 * max(abs(xrange)), max(abs(xrange))) 
    if(is.null(range)) {
      if(min(x) < 0)
        m <- (-1)
      else
        m <- 1
      if(abs(min(x)) > abs(max(x)))
        x <- c(x, abs(min(x)))
      if(abs(max(x)) > abs(min(x)))
        x <- c(x, m * abs(max(x)))
      id <- c(as.character(id), "added")
    } else {
      if(max(range) > max(x)) {
        x <- c(x, max(range))
        id <- c(as.character(id), "added")
      } else x[x > max(range)] <- max(range)
      if(min(range) < min(x)) {
        x <- c(x, min(range))
        id <- c(as.character(id), "added")
      } else x[x < min(range)] <- min(range)
    }
  }

  return(list(id = as.character(id), x = x, range = xrange))
}

