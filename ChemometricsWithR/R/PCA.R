PCA <- function(X, warn = TRUE)
{
  ndf <- nrow(X) - 1
  X.svd <- svd(X)
  varnames <- colnames(X)
  if (is.null(varnames)) varnames <- paste("Var", 1:ncol(X))

  object <- list(scores = X.svd$u %*% diag(X.svd$d),
                 loadings = X.svd$v,
                 var = X.svd$d^2 / ndf,
                 totalvar = sum(X.svd$d^2)/ndf)
  dimnames(object$scores) <- list(rownames(X),
                                  paste("PC", 1:ncol(object$scores)))
  dimnames(object$loadings) <- list(varnames,
                                    paste("PC", 1:ncol(object$loadings)))
  names(object$var) <- paste("PC", 1:length(X.svd$d))

  if (!isTRUE(all.equal(colMeans(X), rep(0, ncol(X)),
                        check.attributes = FALSE))) {
    if (warn)
      warning("Performing PCA on a non-meancentered data matrix!")
    object$centered.data <- FALSE
  } else {
    object$centered.data <- TRUE
  }

  class(object) <- "PCA"

  object
}

variances <- function(object, npc = maxpc)
{
  maxpc <- max(ncol(object$loadings), ncol(object$scores))
  if (npc > maxpc) {
    warning(paste("Maximal number of PCs:", maxpc))
    npc <- maxpc
  }
  
  object$var[1:npc]
}

screeplot <- function(object, type = c("scree", "percentage"), npc, ...)
{
  if (missing(npc)) npc <- length(variances(object))
  type <- match.arg(type)

  vars <- switch(type,
                 scree = log(variances(object)[1:npc]),
                 percentage = cumsum(100*object$var[1:npc]/object$totalvar))
  ylab <- switch(type,
                 scree = "log(variance)",
                 percentage = "% variance")

  barplot(vars, names.arg = 1:npc, xlab = "# PCs", ylab = ylab, ...)
}

reconstruct <- function(object, npc = maxpc)
{
  maxpc <- max(ncol(object$loadings), ncol(object$scores))
  orig <- tcrossprod(scores(object), loadings(object))
  reconstruction <- tcrossprod(scores(object, npc), loadings(object, npc))

  error <- sqrt(mean((orig - reconstruction)^2))
  
  list(reconstruction = reconstruction,
       error = error)
}

project <- function(object, npc = maxpc, newdata, ldngs)  
{
  if (missing(ldngs)) {
    maxpc <- max(ncol(object$loadings), ncol(object$scores))
    ldngs <- loadings(object, npc)
  }
  
  newdata %*% ldngs[, 1:npc]
}

## generic functions... 
scoreplot <- function(object, ...) UseMethod("scoreplot")
loadingplot <- function(object, ...) UseMethod("loadingplot")

scores.PCA <- function(object, npc = maxpc, ...)
{
  maxpc <- max(ncol(object$loadings), ncol(object$scores))
  if (npc > maxpc) {
    warning(paste("Maximal number of PCs:", maxpc))
    npc <- maxpc
  }
  
  object$scores[, 1:npc, drop = FALSE]
}

scoreplot.PCA <- function(object, pc = c(1,2),
                          pcscores = scores(object),
                          show.names = FALSE,
                          xlab, ylab, xlim, ylim, ...)
{
  if (is.null(varnames <- rownames(pcscores)))
    show.names <- FALSE

  ## Preliminaries: determine range of the plot
  if (missing(xlim)) {
    xlim <- unsigned.range(pcscores[,pc[1]])
    if (show.names) xlim <- xlim * 1.2 # more space needed for text
  }
  if (missing(ylim)) {
    ylim <- unsigned.range(pcscores[,pc[2]])
    if (show.names) ylim <- ylim * 1.1
  }

  ## determine default axis labels
  if (!missing(object) && object$centered.data) {
    if (missing(xlab))
      xlab <- paste("PC", pc[1],
                    paste("(",
                          formatC(100 * object$var[pc[1]] / object$totalvar,
                                  format = "f", digits = 1), "%)", sep = ""))
    if (missing(ylab))
      ylab <- paste("PC", pc[2],
                    paste("(",
                          formatC(100 * object$var[pc[2]] / object$totalvar,
                                  format = "f", digits = 1), "%)", sep = ""))
  } else {
    if (missing(xlab)) xlab <- paste("PC", pc[1])
    if (missing(ylab)) ylab <- paste("PC", pc[2])
  }
  
  on.exit(par(op))
  op <- par(pty = "s") # make a square rather than oblong plot

  ## Go!
  plot(pcscores[,pc], type = "n", xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, ...)
  abline(h = 0, v = 0, col = "gray", lty = 2)

  if (show.names) {
    text(pcscores[,pc[1]], pcscores[,pc[2]], varnames, ...)
  } else {
    points(pcscores[,pc], ...)
  }
}
    
loadings.PCA <- function(object, npc = maxpc, ...)
{
  maxpc <- max(ncol(object$loadings), ncol(object$scores))
  if (npc > maxpc) {
    warning(paste("Maximal number of PCs:", maxpc))
    npc <- maxpc
  }
  
  object$loadings[, 1:npc, drop = FALSE]
}

loadingplot.PCA <- function(object, pc = c(1,2),
                            pcloadings = loadings(object),
                            scalefactor = 1,
                            add = FALSE, show.names = FALSE,
                            xlab, ylab, xlim, ylim, col = "blue",
                            min.length = .01, varnames = NULL, ...)
{
  if (add) pcloadings <- pcloadings / scalefactor

  if (is.null(varnames))
    if (is.null(varnames <- rownames(pcloadings)))
      show.names <- FALSE

  ## Preliminaries: determine range of the plot
  if (missing(xlim)) xlim <- unsigned.range(pcloadings[,pc[1]])
  if (missing(ylim)) ylim <- unsigned.range(pcloadings[,pc[2]])

  if (!add) {
    if (show.names) { # may need extra space for names
      xlim <- xlim * 1.2
      ylim <- ylim * 1.1 # no word length to take into account
    }

    ## determine default axis labels
    if (!missing(object) && object$centered.data) {
      if (missing(xlab))
        xlab <- paste("PC", pc[1],
                      paste("(",
                            formatC(100 * object$var[pc[1]] / object$totalvar,
                                    format = "f", digits = 1), "%)", sep = ""))
      if (missing(ylab))
        ylab <- paste("PC", pc[2],
                      paste("(",
                            formatC(100 * object$var[pc[2]] / object$totalvar,
                                    format = "f", digits = 1), "%)", sep = ""))
    } else {
      if (missing(xlab)) xlab <- paste("PC", pc[1])
      if (missing(ylab)) ylab <- paste("PC", pc[2])
    }
     
    on.exit(par(op))
    op <- par(pty = "s") # make a square rather than oblong plot
    
    plot(pcloadings[,pc], type = "n",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }

  nonzeros <- apply(pcloadings[,pc], 1,
                    function(x) sum(x^2) > min.length)
  if (length(col) == 1)
    col <- rep(col, nrow(pcloadings))

  if (show.names) {
    text(pcloadings[nonzeros,pc[1]] * 1.1,
         pcloadings[nonzeros,pc[2]] * 1.1,
         varnames[nonzeros], col = col[nonzeros])
  } else {
    points(pcloadings[,pc], col = col[nonzeros], ...)
  }

  origin <- rep(0, nrow(pcloadings))
  arrows(origin[nonzeros],
         origin[nonzeros],
         pcloadings[nonzeros, pc[1]],
         pcloadings[nonzeros, pc[2]],
         angle = 15, length = .15,
         col = col[nonzeros], ...)

  if (add) { 
    ## biplot: reset the axis system to the loadings plot
    par(new = TRUE)
    plot(pcloadings[,pc], axes = FALSE, type = "n",
         xlim = xlim, ylim = ylim,
         xlab = "", ylab = "")
    axis(3, col = col[1])
    axis(4, col = col[1])
    box(col = 1)
  }
}


biplot.PCA <- function(x, pc = c(1,2),
                       show.names = c("none", "scores", "loadings", "both"),
                       score.col = 1, loading.col = "blue",
                       min.length = .01, varnames = NULL, ...)  
{
  show.names <- match.arg(show.names)
  show.sc.names <- (show.names == "scores" | show.names == "both") 
  show.ld.names <- (show.names == "loadings" | show.names == "both") 
  
  rangx1 <- unsigned.range(x$scores[, pc[1]])
  rangx2 <- unsigned.range(x$scores[, pc[2]])
  rangy1 <- unsigned.range(x$loadings[, pc[1]])
  rangy2 <- unsigned.range(x$loadings[, pc[2]])
  ratio <- max(rangy1/rangx1, rangy2/rangx2)
  xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)

  scoreplot(x, pc = pc, xlim = xlim, ylim = ylim,
            show.names = show.sc.names, col = score.col, ...)
  loadingplot(x, pc = pc, add = TRUE, scalefactor = ratio,
              xlim = xlim*ratio, ylim = ylim*ratio, 
              pch = 2, show.names = show.ld.names, col = loading.col,
              min.length = min.length, varnames = varnames, ...)
}


summary.PCA <- function(object, varperc = 90,
                        pc.select = c(1:5,10), ...)
{
  nr <- nrow(object$scores)
  nc <- nrow(object$loadings)
  npc <- length(object$var)
  
  cat("\nPCA model of a",
      ifelse (object$centered.data, "mean-centered", ""),
      "matrix of", nr, "by", nc)
  
  vartab <- 100 * cbind(object$var, cumsum(object$var)) / object$totalvar
  colnames(vartab) <- c("Var", "Cumul. var.")#)
  npcs <- sum(vartab[,2] < varperc)
  cat("\nNumber of PCs to cover", varperc,
      "percent of the variance:", npcs + 1)
  
  cat("\n\n")
  which.pcs <- pc.select[pc.select <= npc]
  print(vartab[which.pcs,], ...)

  invisible()                        
}
