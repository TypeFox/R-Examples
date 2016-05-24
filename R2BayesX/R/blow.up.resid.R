blow.up.resid <- function(data, x, xnam, response, eta, dimx, cx)
{
  if(!is.null(data) && !is.null(response) && !is.null(eta) && !is.factor(response)) {
    x <- x[order(x[,1L]), ]
    if(!is.matrix(x))
      x <- matrix(x, nrow = 1)
    id <- NULL
    if(cx == "geo.bayesx") {
      co <- x[,2L:3L]
      id <- x[,1L]
    } else co <- matrix(x[,1L:dimx], ncol = dimx)
    xtmp <- data[[xnam[1L]]]
    if(!is.null(xtmp)) {
      ox <- order(xtmp)
      xtmp <- xtmp[ox]
      response <- response[ox]
      if(!is.matrix(eta))
        eta <- matrix(eta, ncol = 1L)
      eta <- eta[ox, ]
      ind <- unique.id(xtmp)
      if(max(ind) <= nrow(x)) {
        x <- as.data.frame(x[ind,])
        x$pcat80 <- NULL
        x$pcat95 <- NULL
        x$pcat80_sim <- NULL
        x$pcat95_sim <- NULL
        x$pcat80_tot <- NULL
        x$pcat95_tot <- NULL
        x$pcat80tot_sim <- NULL
        x$pcat95tot_sim <- NULL
        for(k in 1L:length(xnam))
          eval(parse(text = paste("x$\'", xnam[k] , "\' <- NULL", sep = "")))
        x <- as.matrix(x)
        if(!is.matrix(eta))
          eta <- matrix(eta, ncol = 1L)
        if(length(xm <- colnames(x)[colnames(x) == "pmean"]) < 1L)
          xm <- 1L
        pres <- response - eta[,1L] + x[,xm]
        x <- cbind(co[ind, ], pres, id[ind])
        if(ncol(x) < 3L)
          colnames(x) <- c("x.co", "partial.resids")
        else {
          if(!is.null(id))
            colnames(x) <- c("x.co", "y.co", "partial.resids", "id")
          else
            colnames(x) <- c("x.co", "y.co", "partial.resids")
        }
        rownames(x) <- 1L:nrow(x)
      }
    }
  }

  return(x)
}

