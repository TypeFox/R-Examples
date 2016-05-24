temp_capture.output <- function (..., file = NULL, append = FALSE, type = c("output",
                                                             "message"), split = FALSE) {
  args <- substitute(list(...))[-1L]
  type <- match.arg(type)
  rval <- NULL
  closeit <- TRUE
  if (is.null(file))
    file <- textConnection("rval", "w", local = TRUE)
  else if (is.character(file))
    file <- file(file, if (append)
      "a"
      else "w")
  else if (inherits(file, "connection")) {
    if (!isOpen(file))
      open(file, if (append)
        "a"
        else "w")
    else closeit <- FALSE
  }
  else stop("'file' must be NULL, a character string or a connection")
  sink(file, type = type, split = split)
  on.exit({
    sink(type = type, split = split)
    if (closeit) close(file)
  })
  pf <- parent.frame()
  evalVis <- function(expr) withVisible(eval(expr, pf))
  for (i in seq_along(args)) {
    expr <- args[[i]]
    tmp <- switch(mode(expr), expression = lapply(expr, evalVis),
                  call = , name = list(evalVis(expr)), stop("bad argument"))
    for (item in tmp) if (item$visible)
      print(item$value)
  }
  on.exit()
  sink(type = type, split = split)
  if (closeit)
    close(file)
  if (is.null(rval))
    invisible(NULL)
  else rval
}

volTriangulationWrapper <- function(vertices) {
  if (sessionInfo()$R.version$`svn rev` < "69993") {
    ## FR->FR should be obsolete some day... load def of capture.output from R devel, future 3.3.0
    capture.output <- temp_capture.output
  } ## else R already has the right capture.output
  abyss <- capture.output(vT <- try(volTriangulation(vertices)),type="message")
  if (inherits(vT,"try-error")) {
    ## **FR->FR** here because we don't use rcdd in spaMM
    cnames <- colnames(vertices) ## redundant loses names
    vertices <- q2d(redundant(d2q(cbind(0, 1, as.matrix(vertices))), representation="V")$output[, -c(1:2), drop=FALSE]) ## FR->FR heavy solution
    if (nrow(vertices)<=ncol(vertices)) {
      return(try(stop("nrow <= ncol in volTriangulationWrapper() for the vertices\n of the convex hull of the input points.",call.=FALSE),silent=TRUE))
    }
    origin <- vertices[1,]
    DV <- sweep(vertices[-1,,drop=FALSE], 2, origin, `-`)
    orthog <- qr.Q(qr(t(DV))) ## orthonormal basis of dim = nrow(vertices) possibly lower than ncol(vertices)
    projcoefs <- rbind(0,DV %*% orthog) ## add coeffs of projected origin
    vT <- volTriangulation(projcoefs) ## triangulation on projected coordinates
    ## back to original coordinates (volumes unchanged by orthog projection)
    vT$bary <- sweep(vT$bary %*% t(orthog), 2, origin, `+`)
    vT$vertices <- sweep(vT$vertices %*% t(orthog), 2, origin, `+`)
    colnames(vT$vertices) <- cnames
  }
  return(vT) ## vT$vertices should always have colnames
}
