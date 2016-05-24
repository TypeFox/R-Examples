likertMosaic <- function(x, ...)
  UseMethod("likertMosaic")

likertMosaic.formula <- function(x, data, ReferenceZero=NULL, spacing=NULL,
                                 ..., between.y=c(1.2, .3)) {
  spacing.in <- spacing
  varNamesUsed <- getVarNames(x, data)
  ## list(QuestionName, CondNames, LevelNames) ## Subset of data columns actually used

  data.list <- getLikertData(data, varNamesUsed) ## list(Question, Conditions, Nums)
  Nums <- data.list$Nums
  rownames(Nums) <- data.list$Question[[1]]

  if (!is.null(varNamesUsed$CondNames)) {
    groupLengths <- table(data.list$Conditions[[1]])
    numberofbetweens <- as.vector(rbind(1, groupLengths-1))
    between.index <- rep(rep(c(1,2), length=length(numberofbetweens)), numberofbetweens)[-1]
    colorset <- ColorSet(ncol(Nums), ReferenceZero)
    n.likert.levels <- length(colorset) + (0 %in% colorset) + 2 ## extra two for left and right padding
    spacing=c(list(unit(between.y[between.index], "lines")),
      spacing_highlighting()(n.likert.levels))
  } else
    spacing=spacing_highlighting

  if (!is.null(spacing.in)) spacing <- spacing.in

  likertMosaic(Nums, ReferenceZero=ReferenceZero, ...,
               spacing=spacing,
               Conditions=data.list$Conditions)
}

likertMosaic.default <- function(x, ...) { ## most likely for a vector
  x.input <- x
  if (is.null(dim(x))) {
    x <- t(x)
    if (is.null(dimnames(x))) dimnames(x) <- list("", letters[seq(along=x)])
    dimnames(x)[[1]] <- ""
  }
  likertMosaic(x)
}

likertMosaic.list <- function(x, ...) {
  lapply(x, likertMosaic, ...)
}


likertMosaic.data.frame <- function(x, ...) {
  likertMosaic(data.matrix(x), ...)
}


likertMosaic.matrix <- function(x, ...,
                                split_vertical=c(FALSE,TRUE),
                                rot_labels=c(90,0,0,0),
                                just_labels=c("left","center","center","right"),
                                labels=c(TRUE,FALSE)) {
  if (is.null(dimnames(x)))
    dimnames(x) <- list(NULL, letters[seq(length=ncol(x))])
  if (is.null(dimnames(x)[[2]]))
    dimnames(x)[[2]] <- letters[seq(length=ncol(x))]
  likertMosaic.array(x,
                     split_vertical=split_vertical,
                     rot_labels=rot_labels,
                     just_labels=just_labels,
                     labels=labels,
                     ...)
}

likertMosaic.array <- function(x, ReferenceZero=NULL, col=NULL, main=NULL, ...,
                               as.percent=FALSE,
                               variable.width=NULL,
                               positive.order=FALSE,
                               Conditions=NULL,
                               x.legend=list(text=list(dimnames(x)[[ndim]]), columns=x.dim[ndim],
                                 space="bottom", size=2, cex=.8, between=.6,
                                 rect=list(col=col, border="white")),
                               legend.y=0.05,
                               ## arguments following are mosaic or strucplot arguments
                               spacing=spacing_highlighting,
                               split_vertical=c(TRUE,FALSE),
                               margins=c(3,2,4,22), ## clockwise from top
                               keep_aspect=FALSE,
                               rot_labels=c(0,0,90,0),
                               just_labels=c("center","center","center","right"),
                               labels=c(TRUE,TRUE,FALSE,TRUE),
                               varnames=FALSE,
                               zero_size=0,
                               gp=gpar(fill=col.extended, ## fill color for tiles
                                 col=0),                   ## border color for tiles
                               colorFunction="diverge_hcl",
                               colorFunctionOption="lighter"
                               ) {

  x.dimnames <- dimnames(x)
  x.dim <- dim(x)
  ndim <- length(x.dimnames)

  xmat <- x
  if (ndim > 2) {
    dim(xmat) <- c(prod(x.dim[-ndim]), x.dim[ndim])
    dimnames(xmat) <- list(NULL, x.dimnames[[ndim]])
  }

  if (as.percent != FALSE) {
    rsx <- rowSums(xmat)
    xmat <- xmat / rsx
  }

  xmat.lik <- as.likert(xmat, ReferenceZero=ReferenceZero, padding=TRUE, reverse.left=FALSE)
  attr.xmat.lik <- attributes(xmat.lik)

  if (as.percent != FALSE && !is.null(variable.width)) {
    xmat.lik <- xmat.lik * rsx
  }

  if (positive.order) {
    if (is.null(Conditions) || ncol(Conditions)==0)
      xmat.lik <- xmat.lik[rev(attr.xmat.lik$positive.order), , drop=FALSE]
    else
      xmat.lik <- xmat.lik[order(Conditions$Subtable, order(rev(attr.xmat.lik$positive.order))), , drop=FALSE]
  }

  if (ndim > 2) {
    xmat.lik.names.3 <- dimnames(xmat.lik)[[2]]
    dim(xmat.lik) <- c(x.dim[-ndim], length(xmat.lik.names.3))
    dimnames(xmat.lik) <- c(x.dimnames[-ndim], Levels=list(xmat.lik.names.3))
  }

  ## vcd::mosaic(xmat.lik, split_vertical=c(TRUE,FALSE))  ## the winner (for 2)
  ## vcd::mosaic(xmat.lik, split_vertical=c(TRUE,FALSE,TRUE))  ## the winner (for 3)

  if (is.null(col))
    col <- likertColor(attr.xmat.lik$nlevels,
                       ReferenceZero=ReferenceZero,
                       colorFunction=colorFunction,
                       colorFunctionOption=colorFunctionOption)
  col.extended=c(
    "transparent",
    col[attr.xmat.lik$color.seq],
    "transparent")
  dim(col.extended) <- c(1, length(col.extended)) ## this line needed with vcd_1.2-13 and earlier

  result <-
    vcd::mosaic(xmat.lik,
                keep_aspect=keep_aspect,
                spacing=spacing,
                split_vertical=split_vertical,
                rot_labels=rot_labels,
                just_labels=just_labels,
                labels=labels,
                varnames=varnames,
                margins=margins, ## clockwise from top
                main=main,
                zero_size=0,
                gp=gp,
                ...)

  lattice::draw.key(x.legend, draw=TRUE, vp=viewport(x=.5, y=legend.y))

  invisible(result)
}
