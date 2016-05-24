venneuler <- function(combinations, weights, ...) {
  if (missing(combinations)) stop("combinations must be specified")
  if (inherits(combinations, "table")) {
    if (!missing(weights)) warning("combinations is a table yet weights are also specified - ignoring weights")
    weights <- as.vector(combinations)
    rnm <- rep(rownames(combinations), dim(t)[2])
    cnm <- rep(colnames(combinations), each=dim(t)[1])
    names(weights) <- paste(rnm, cnm, sep="&")
    if (all(weights == 0)) stop("all weights are zero")
    weights <- weights[weights != 0]
    combinations <- names(weights)
  } else if (missing(weights) && is.numeric(combinations) && is.null(dim(combinations))) {
    if (is.null(names(combinations))) stop("combinations are a numeric vector but without names")
    weights <- combinations
    combinations <- names(combinations)
  }
  if (is.data.frame(combinations)) combinations <- as.matrix(combinations)
  if (is.matrix(combinations) && (is.numeric(combinations) || is.logical(combinations))) {
    if (is.null(colnames(combinations))) colnames(combinations) <- LETTERS[seq.int(dim(combinations)[2])]
    ## aggregate all entries using a hased environment -- we could probably devise a smarter way if we cared ...
    e <- new.env(TRUE, emptyenv())
    cn <- colnames(combinations)
    if (is.logical(combinations)) { for (i in seq.int(dim(combinations)[1])) if (any(combinations[i,])) {
      ec <- paste(cn[combinations[i,]], collapse='&')
      e[[ec]] <- if (is.null(e[[ec]])) 1L else (e[[ec]] + 1L)
    } } else if (is.numeric(combinations)) for (i in seq.int(dim(combinations)[1])) if (any(combinations[i,] != 0)) {
      ec <- paste(cn[combinations[i,] != 0], collapse='&')
      e[[ec]] <- (if (is.null(e[[ec]])) 0 else e[[ec]]) + sum(combinations[i,])
    }
    en <- ls(e, all.names=TRUE)
    weights <- as.numeric(unlist(lapply(en, get, e)))
    combinations <- as.character(en)
  }
  if (is.matrix(combinations) && is.character(combinations) && dim(combinations)[2] == 2) {
    vd <- .jnew("edu/uic/ncdm/venn/data/VennData", as.character(combinations[,1]), as.character(combinations[,2]))
  } else {
    if (!is.character(combinations)) stop("combinations must be either a character vector, a table, a named numeric vector or a character matrix with two columns")
    if (missing(weights)) weights <- rep(1, length(combinations))
    vd <- .jnew("edu/uic/ncdm/venn/data/VennData", as.character(combinations), as.double(weights))
  }
  a <- .jnew("edu/uic/ncdm/venn/VennAnalytic")
  g <- .jcall(a, "Ledu/uic/ncdm/venn/VennDiagram;", "compute", vd)
  ct <- lapply(.jevalArray(.jfield(g, "[[D", "centers", convert=FALSE)), .jevalArray)
  n <- length(ct)
  ct <- matrix(unlist(ct), ncol=2, byrow=TRUE)
  colnames(ct) <- c("x", "y")
  diam <- .jevalArray(.jfield(g, "[D", "diameters", convert=FALSE))
  areas <- .jevalArray(.jfield(g, "[D", "areas", convert=FALSE))
  res <- .jevalArray(.jfield(g, "[D", "residuals", convert=FALSE))
  col <- .jevalArray(.jfield(g, "[D", "colors", convert=FALSE))
  lab <- .jevalArray(.jfield(g, "[Ljava/lang/String;", "circleLabels", convert=FALSE))
  rownames(ct) <- lab
  names(diam) <- lab
  names(col) <- lab
  names(res) <- .jevalArray(.jfield(g, "[Ljava/lang/String;", "residualLabels", convert=FALSE))
  structure(list(centers=ct, diameters=diam, colors=col, labels=lab, residuals=res,
            stress=.jfield(g, "D", "stress"), stress01=.jfield(g, "D", "stress01"),
	    stress05=.jfield(g, "D", "stress05")), class="VennDiagram")
}

## Note: in col.fn we need more croma and less luminance than usual, because we'll be plotting with reduced alpha
plot.VennDiagram <- function(x, col, col.fn = function(col) hcl(col * 360, 130, 60), alpha=0.3, main=NULL, edges=200, border=NA, col.txt=1, ...) {
  # calculate total extents
  xtp <- x$centers + x$diameters / 2
  xtm <- x$centers - x$diameters / 2
  xr <- range(c(xtp[,1], xtm[,1]))
  yr <- range(c(xtp[,2], xtm[,2]))
  # create canvas
  plot.new()
  plot.window(xr, yr, "", asp = 1)
  # adjust alpha for all colors if specified
  n <- length(x$diameters)
  if (missing(col)) col <- col.fn(x$colors)
  if (length(col) < n) col <- rep(col, length.out=n)
  if (!is.na(alpha)) {
    col <- col2rgb(col) / 255
    col <- rgb(col[1,], col[2,], col[3,], alpha)
  }
  # prepare circle coordinates
  s <- seq.int(edges) / edges * 2 * pi
  sx <- cos(s) / 2 # VD uses diameter, not radius
  sy <- sin(s) / 2
  if (!is.null(border)) border <- rep(border, length.out=n)
  # plot all circles
  for (i in seq.int(n))
    polygon(x$centers[i, 1] +  x$diameters[i] * sx, x$centers[i, 2] + x$diameters[i] * sy, col = col[i], border = border[i])
  # if col.txt is not NA, plot the circle text
  if (!all(is.na(col.txt))) text(x$centers, x$labels, col=col.txt)
  # finish with title
  title(main = main, ...)
  invisible(NULL)
}
