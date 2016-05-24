## Assign new data using wcc with unit weight vector

wccassign <- function(x, data)
{
  nobj <- nrow(data)
  nvar <- ncol(data)
  nunits <- nrow(x$grid$pts)
  tw <- x$trwdth
  classif <- rep(0, nobj)
  wccs <- rep(0, nobj)

  if (tw > 0)
    wghts <- 1 - (0:tw)/tw
  else
    wghts <- 1
  
  codes <- t(x$codes)
  acors <- x$acors

  if (missing(data) & !is.null(x$data)) {
    data <- t(x$data)
    data.acors <- x$data.acors
  } else {
    data <- t(data)
    data.acors <- wacmat(data, tw, wghts, do.transpose = FALSE)
  }
  
  res <- .C("wccassign",
            data = as.double(data),
            data.acors = as.double(data.acors),
            codes = as.double(codes),
            acors = as.double(acors),
            as.integer(tw),
            wghts = as.double(wghts),
            classif = as.integer(classif),
            wccs = as.double(wccs),
            as.integer(nobj),
            as.integer(nvar),
            as.integer(nunits),
            PACKAGE = "wccsom")

  list(classif = res$classif, wccs = res$wccs)
}
