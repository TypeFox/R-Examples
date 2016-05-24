density.plot.lambdaMD <- function(x, by = "column", col, 
                                  xlim = c(0,1), ylim,  
                                  main = "", sub = NULL, xlab,
                                  ylab, lty = par("lty"), lwd = par("lwd"), ...) {
  readpars <- par(no.readonly = TRUE)
  if (all(class(x) != "lambdaMD"))
    stop("works only with output from `lambda.MD'")

  getY <- function(x) x[[1]]$y
  get2 <- function(x) x[2]
  
  tnames <- strsplit(colnames(x), "lambda.")
  idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
  idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE,
                                      nrow = length(idx), ncol = 
                                      length(idx[[1]]))))
  idx <- lapply(idx, as.character)
  idx <- lapply(idx, unique)
  lidx <- sapply(idx, length)
  names(lidx) <- names(idx) <- c("rows", "columns")

  if (is.mcmc(x)) { 
    x <- array(t(x), dim = c(sapply(idx, length), nrow(x)),
               dimnames = list(rows = idx[[1]], columns = idx[[2]],
                 simulations = 1:nrow(x)))
  }
  dens <- apply(x, c(1,2), density)
  if (missing(ylim)) { 
    ylim <- apply(apply(dens, c(1,2), getY), c(2,3), max)
  }
  
  if (by == "row") { 
    par(mfrow = c(lidx[1], 1), ...)
    if (missing(xlab)) {
      xlab <- paste("Percentage", unlist(idx[[1]]))
    }
    else if (length(xlab) != length(idx[[1]])) {
      warning(paste("xlab needs to be", lidx[1], "long"))
      xlab <- rep(xlab, lidx[1])[1:lidx[1]]
    }
    if (missing(col)) {
      col <- rainbow(lidx[2])
    }
    else {
      if (length(col) < lidx[2])
        col <- rep(col, idx[[2]])[1:lidx[2]]
      if (length(col) > lidx[2])
        col <- col[1:lidx[2]]
    }
    names(col) <- idx[[2]]
    names(xlab) <- idx[[1]]
    for (ii in idx[[1]]) {
      plot(dens[ii, 1][[1]], type = "l", xlim = xlim, col = col[1],
           ylim = c(0, max(ylim[ii,])), main = main, sub = sub,
           xlab = xlab[ii], ylab = "Density", lty = lty,
           lwd = lwd)
      
      for (jj in idx[[2]][2:lidx[2]])
        lines(dens[ii,jj][[1]], lty = lty, lwd = lwd, col = col[jj])
    }
  }
  if (by == "column") {
    par(mfrow = c(lidx[2], 1), ...)
    if (missing(xlab)) {
      xlab <- paste("Percentage", unlist(idx[[2]]))
    }
    else {
      if (length(xlab) != lidx[1]) {
        warning(paste("xlab needs to be", lidx[1], "long"))
        xlab <- rep(xlab, lidx[2])[1:lidx[2]]
      }
    }
    if (missing(col)) {
      col <- rainbow(lidx[1])
    }
    else {
      if (length(col) < lidx[1])
        col <- rep(col, idx[[1]])[1:lidx[1]]
      if (length(col) > lidx[1])
        col <- col[1:lidx[1]]
    }
    names(col) <- idx[[1]]
    names(xlab) <- idx[[2]]
    for (ii in idx[[2]]) {
      plot(dens[1, ii][[1]], type = "l", xlim = xlim, col = col[1],
           ylim = c(0, max(ylim[,ii])), main = main, sub = sub,
           xlab = xlab[ii], ylab = "Density", lty = lty,
           lwd = lwd)
      for (jj in idx[[1]][2:lidx[1]])
        lines(dens[jj,ii][[1]], lty = lty, lwd = lwd, col = col[jj])
    }
  }
  par(readpars)
  return(invisible(x))
}
