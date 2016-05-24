
# Convert graphics-system-neutral "picture" into grid gList

# Viewport from picture
pictureVP <- function(picture, exp=0.05, xscale=NULL, yscale=NULL,
                      distort=FALSE, ...) {
    if (is.null(xscale) || is.null(yscale)) {
        xscale <- picture@summary@xscale
	yscale <- picture@summary@yscale
    }
    xscale <- xscale + exp*c(-1, 1)*diff(range(xscale))
    yscale <- yscale + exp*c(-1, 1)*diff(range(yscale))
    # If distort=TRUE, having the two layers of viewports is
    # massively redundant, BUT I'm keeping it so that either
    # way there is the same viewport structure, which I think
    # is beneficial if anyone ever wants to make use of
    # these viewports (otherwise they would need to figure
    # out whether a picture grob has one or two viewports).
    vpStack(viewport(name="picture.shape", ...,
                     layout=grid.layout(1, 1,
                       widths=abs(diff(xscale)),
                       heights=abs(diff(yscale)),
                       respect=!distort)),
            viewport(name="picture.scale",
                     layout.pos.col=1,
                     xscale=xscale,
                     yscale=yscale))
}

##################
# Explode a PictureOp whose path has more than one "move"
# into multiple PictureOps
##################
explodePath <- function(path, fill) {
    ops <- attr(path@x, "names")
    if (length(path@x) > 1) {
        if (length(ops) > 1) {
            moves <- grep("move", ops)
            npaths <- length(moves)
            if (npaths > 1) {
                newpaths <- vector("list", npaths)
                for (i in 1:(npaths - 1)) {
                    index <- moves[i]:(moves[i + 1] - 1)
                    if (length(index) > 1) {
                        if (fill) {
                            newpaths[[i]] <- new("PictureFill",
                                                 x=path@x[index],
                                                 y=path@y[index],
                                                 rule="winding",
                                                 lwd=path@lwd,
                                                 lty=path@lty,
                                                 rgb=path@rgb,
                                                 lineend=path@lineend,
                                                 linejoin=path@linejoin,
                                                 linemitre=path@linemitre)
                        } else {
                            newpaths[[i]] <- new("PictureStroke",
                                                 x=path@x[index],
                                                 y=path@y[index],
                                                 lwd=path@lwd,
                                                 lty=path@lty,
                                                 rgb=path@rgb,
                                                 lineend=path@lineend,
                                                 linejoin=path@linejoin,
                                                 linemitre=path@linemitre)
                        }
                    }
                }
                index <- moves[npaths]:length(ops)
                if (length(index) > 1) {
                    if (fill) {
                        newpaths[[npaths]] <-
                            new("PictureFill",
                                x=path@x[moves[npaths]:length(ops)],
                                y=path@y[moves[npaths]:length(ops)],
                                rule="winding",
                                lwd=path@lwd,
                                lty=path@lty,
                                rgb=path@rgb,
                                lineend=path@lineend,
                                linejoin=path@linejoin,
                                linemitre=path@linemitre)
                    } else {
                        newpaths[[npaths]] <- 
                            new("PictureStroke",
                                x=path@x[moves[npaths]:length(ops)],
                                y=path@y[moves[npaths]:length(ops)],
                                lwd=path@lwd,
                                lty=path@lty,
                                rgb=path@rgb,
                                lineend=path@lineend,
                                linejoin=path@linejoin,
                                linemitre=path@linemitre)
                    }
                }
                newpaths[!sapply(newpaths, is.null)]
            } else {
                path
            }
        } else {
            # If there are 'x' values BUT no 'ops' just return the path
            path
        }
    } else {
        # If there are no 'x' values, return an empty list
        list()
    }
}

setGeneric("explode",
           function(object, ...) {
               standardGeneric("explode")
           })

setMethod("explode", signature(object="PictureStroke"),
          function(object, ...) {
              explodePath(object, FALSE)
          })

setMethod("explode", signature(object="PictureFill"),
          function(object, ...) {
              explodePath(object, TRUE)
          })

# Exploding a PictureText does nothing to the object
setMethod("explode", signature(object="PictureText"),
          function(object, ...) {
              object
          })

##################
# Modify the paths in a PictureChar (known to be a traced character):
# First explode the character into multiple paths.
# If stroking, add path start at end to close each path.
# If filling, fill first path with path colour, but all
# subsequent paths with a background colour
# (simple heuristic that may or may not work).
##################
fixPath <- function(path, i, fill, bg) {
    # At this point I know that each path has only one "move"
    # AND that each path has at least two locations
    # (these are sorted out by explodePath()
    if (fill) {
        if (i == 1) {
            path
        } else {
            new("PictureFill",
                x=path@x,
                y=path@y,
                lwd=path@lwd, lty=path@lty, rgb=bg,
                lineend=path@lineend, linejoin=path@linejoin,
                linemitre=path@linemitre)
        }
    } else {
        new("PictureStroke",
            x=c(path@x, path@x[1]),
            y=c(path@y, path@y[1]),
            # Use a light stroke
            lwd=path@lwd, lty=path@lty, rgb=path@rgb,
            lineend=path@lineend, linejoin=path@linejoin,
            linemitre=path@linemitre)
    }
}

setMethod("explode", signature(object="PictureChar"),
          function(object, fill, bg, ...) {
              paths <- explodePath(object, fill)
              np <- length(paths)
              if (np > 0) {
                  # bg can be a named vector
                  bgLetters <- names(bg)
                  if (length(bgLetters) > 0) {
                      letterBG <- bg[object@char]
                      if (!is.na(letterBG))
                          bg <- letterBG
                  }
                  newpaths <- vector("list", np)
                  for (i in 1:np) {
                      newpaths[[i]] <- fixPath(paths[[i]], i, fill, bg)
                  }
                  paths <- newpaths
              }
              paths
          })

.bgText.default <- c(a="white",
                     b="white",
                     d="white",
                     e="white",
                     g="white",
                     i="black",
                     j="black",
                     o="white",
                     p="white",
                     q="white",
                     A="white",
                     B="white",
                     D="white",
                     O="white",
                     P="white",
                     Q="white",
                     R="white")


##################
# Convert picture or path into single grob
# For using picture as a one-off (e.g., plot background)
##################
# Generic grobify() function
setGeneric("grobify",
           function(object, ...) {
               standardGeneric("grobify")
           })

picStrokeGrob <- function(...) {
    gTree(..., cl="picstroke")
}

makeContent.picstroke <- function(x) {
    # Figure out what lwd and lty really are
    lwd <- convertWidth(unit(x$lwd, "native"), "bigpts", valueOnly=TRUE)
    lty <- fixLTY(x$lty, x$lwd)
    child <- polylineGrob(x$x, x$y,
                          default.units=x$default.units,
                          id.lengths=x$id.lengths,
                          gp=gpar(lwd=lwd, lty=lty, col=x$col, fill=NA,
                                  lineend=x$lineend, linejoin=x$linejoin,
                                  linemitre=x$linemitre))
    setChildren(x, gList(child))
}

# Individual path converted into grob
setMethod("grobify", signature(object="PictureStroke"),
          function(object, ..., fillText, bgText, sizeByWidth, use.gc=TRUE) {
              if (length(object@x) > 1) {
                  paths <- explode(object)
                  if (is.list(paths)) {
                      pathX <- lapply(paths, slot, "x")
                      pathY <- lapply(paths, slot, "y")
                  } else {
                      pathX <- list(paths@x)
                      pathY <- list(paths@y)
                  }
                  if (use.gc) {
                      picStrokeGrob(x=unlist(pathX),
                                    y=unlist(pathY),
                                    default.units="native",
                                    id.lengths=sapply(pathX, length),
                                    lwd=object@lwd,
                                    lty=object@lty,
                                    col=object@rgb,
                                    lineend=object@lineend,
                                    linejoin=object@linejoin,
                                    linemitre=object@linemitre,
                                    ...)
                  } else {
                      polylineGrob(x=unlist(pathX),
                                   y=unlist(pathY),
                                   default.units="native",
                                   id.lengths=sapply(pathX, length),
                                   ...)
                  }
              } else {
                  NULL
              }
          })

setMethod("grobify", signature(object="PictureFill"),
          function(object, ..., fillText, bgText, sizeByWidth, use.gc=TRUE) {
              if (length(object@x) > 1) {
                  paths <- explode(object)
                  if (is.list(paths)) {
                      pathX <- lapply(paths, slot, "x")
                      pathY <- lapply(paths, slot, "y")
                  } else {
                      pathX <- list(paths@x)
                      pathY <- list(paths@y)
                  }
                  if (use.gc) {
                      pathGrob(x=unlist(pathX),
                               y=unlist(pathY),
                               default.units="native",
                               id.lengths=sapply(pathX, length),
                               rule=switch(object@rule,
                                 nonzero="winding", "evenodd"),
                               gp=gpar(col=NA, fill=object@rgb,
                                       lineend=object@lineend,
                                       linejoin=object@linejoin,
                                       linemitre=object@linemitre),
                               ...)
                  } else {
                      pathGrob(x=unlist(pathX),
                               y=unlist(pathY),
                               default.units="native",
                               id.lengths=sapply(pathX, length),
                               rule=switch(object@rule,
                                 nonzero="winding", "evenodd"),
                               ...)
                  }
              } else {
                  NULL
              }
          })

setMethod("grobify", signature(object="PictureText"),
          function(object, FUN=grobify, ..., 
                   fillText=FALSE, bgText=.bgText.default,
                   sizeByWidth=TRUE, use.gc=TRUE) {
              if (length(object@letters) == 0) {
                  if (use.gc) {
                      pictureTextGrob(object@string,
                                      object@x, object@y,
                                      object@w, object@h,
                                      object@angle,
                                      object@letters,
                                      ...,
                                      sizeByWidth=sizeByWidth,
                                      gp=gpar(col=object@rgb))
                  } else {
                      pictureTextGrob(object@string,
                                      object@x, object@y,
                                      object@w, object@h,
                                      object@angle,
                                      object@letters,
                                      ...,
                                      sizeByWidth=sizeByWidth)
                  }
              } else {
                  gTree(string=as.character(object@string), 
                        children=do.call("gList",
                          lapply(object@letters, FUN=FUN, ...,
                                 fillText=fillText, bgText=bgText,
                                 sizeByWidth=sizeByWidth, use.gc=use.gc)),
                        cl="pictureletters")
              }
          })

# PictureChars no longer need to be exploded and drawn as
# separate paths because can now draw them as single path (with holes)
setMethod("grobify", signature(object="PictureChar"),
          function(object, ...,
                   fillText=FALSE, bgText=NA,
                   sizeByWidth=TRUE, use.gc=TRUE) {
              paths <- explode(object, FALSE, NA)
              if (length(paths) > 0) {
                  pathX <- lapply(paths, slot, "x")
                  pathY <- lapply(paths, slot, "y")
                  if (use.gc) {
                      pathGrob(x=unlist(pathX),
                               y=unlist(pathY),
                               default.units="native",
                               id.lengths=sapply(pathX, length),
                               rule="winding",
                               gp=gpar(col=NA, fill=object@rgb,
                                       lineend=object@lineend,
                                       linejoin=object@linejoin,
                                       linemitre=object@linemitre),
                               ...)
                  } else {
                      pathGrob(x=unlist(pathX),
                               y=unlist(pathY),
                               default.units="native",
                               id.lengths=sapply(pathX, length),
                               rule="winding",
                               ...)
                  }
              } else {
                  NULL
              }
          })

pictureHull <- function(object) {
    allx <- unlist(lapply(object@paths, function(x) { x@x }))
    ally <- unlist(lapply(object@paths, function(x) { x@y }))
    hull <- chull(allx, ally)
    list(x=allx[hull], y=ally[hull])
}

setMethod("grobify", signature(object="Picture"),
          function(object,
                   x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                   width=unit(1, "npc"), height=unit(1, "npc"),
                   just="centre", xscale=NULL, yscale=NULL, exp=0.05,
                   FUN=grobify, distort=FALSE, ..., name=name, gp=gpar()) {
              gTree(childrenvp=pictureVP(object,
                      exp=exp, xscale=xscale, yscale=yscale,
                      distort=distort,
                      x=x, y=y, width=width, height=height,
                      just=just, gp=gp),
                    children=do.call("gList",
                      lapply(object@paths, FUN=FUN,
                             vp=vpPath("picture.shape", "picture.scale"),
                             ...)),
                    # FIXME: editDetails method for "picture" class
                    # to update these bounds if picture grob edited?
                    # (or at least warn that bounds are not updated)
                    hull=pictureHull(object),
                    name=name, cl="picture")
          })

# grobX, grobY methods for gTree of class "picture"
xDetails.picture <- function(x, theta) {
    # Generate a polygon based on picture convex hull
    # NOTE: make sure it has same vpPath as actual paths
    grobX(polygonGrob(x$hull$x, x$hull$y,
                      default.units="native",
                      vp=vpPath("picture.shape", "picture.scale")),
          theta)
}

yDetails.picture <- function(x, theta) {
    # Generate a polygon based on picture convex hull
    # NOTE: make sure it has same vpPath as actual paths
    grobY(polygonGrob(x$hull$x, x$hull$y,
                      default.units="native",
                      vp=vpPath("picture.shape", "picture.scale")),
          theta)
}

##################
# Convert picture or path into "multiple" grob with no viewports
# For using picture as data symbol
##################
setGeneric("symbolize",
           function(object,
                    x=unit(0.5, "npc"),
                    y=unit(0.5, "npc"),
                    size=unit(1, "npc"),
                    units="npc", ...) {
               standardGeneric("symbolize")
           })

symbolLocn <- function(object, x, y, size, units,
                       xscale, yscale) {
    n <- max(length(x), length(y))
    x <- rep(x, length.out=n)
    y <- rep(y, length.out=n)
    size <- rep(size, length.out=n)
    if (!is.unit(x))
        x <- unit(x, units)
    if (!is.unit(y))
        y <- unit(y, units)
    if (!is.unit(size))
        size <- unit(size, units)
    # Do everything in INCHES to avoid too much calculation on units
    x <- convertX(x, "inches", valueOnly=TRUE)
    y <- convertY(y, "inches", valueOnly=TRUE)
    sizew <- convertWidth(size, "inches", valueOnly=TRUE)
    sizeh <- convertHeight(size, "inches", valueOnly=TRUE)
    if (is.null(xscale))
        rx <- range(object@x)
    else
        rx <- xscale
    if (is.null(yscale))
        ry <- range(object@y)
    else
        ry <- yscale
    # Determine width and height from size and scale
    # so that the largest scale dimension
    # is given the smallest physical dimension
    sizeAspect <- sizeh / sizew
    scaleAspect <- diff(ry) / diff(rx)
    width <- ifelse(scaleAspect < sizeAspect,
                    sizew, sizeh/scaleAspect)
    height <- ifelse(scaleAspect < sizeAspect,
                     sizew*scaleAspect, sizeh)
    lwd <- width*object@lwd/diff(rx)*72
    # Scale object@x/y to [-0.5, 0.5]
    # and then multiply by width/height
    wx <- rep((object@x - mean(rx))/abs(diff(rx)), n)*
        rep(width, each=length(object@x))
    hy <- rep((object@y - mean(ry))/abs(diff(ry)), n)*
        rep(height, each=length(object@y))
    # Replicate x/y by length object@x/y
    # NOTE object@x and object@y have same length
    xx <- rep(x, rep(length(object@x), n))
    yy <- rep(y, rep(length(object@x), n))
    list(x=xx + wx, y=yy + hy, n=n, lwd=lwd)
}

makeContent.symbolStroke <- function(x) {
    locn <- symbolLocn(x$object, x$x, x$y, x$size,
                       x$units, x$xscale, x$yscale)
    # Create id to distinguish separate symbols
    id <- rep(1:locn$n, each=length(x$object@x))
    # Generate grob representing symbols
    if (x$use.gc) {
        lwd <- locn$lwd
        lty <- fixLTY(x$object@lty, x$object@lwd)
        child <- do.call("polylineGrob",
                         c(list(x=locn$x, y=locn$y, id=id,
                                default.units="inches",
                                gp=gpar(lwd=lwd,
                                    lty=lty,
                                    col=x$object@rgb,
                                    lineend=x$object@lineend,
                                    linejoin=x$object@linejoin,
                                    linemitre=x$object@linemitre)),
                           x$poly.args))
    } else {
        child <- do.call("polylineGrob",
                         c(list(x=locn$x, y=locn$y, id=id,
                                default.units="inches", vp=x$vp),
                           x$poly.args))
    }
    setChildren(x, gList(child))
}

makeContent.symbolFill <- function(x) {
    locn <- symbolLocn(x$object, x$x, x$y, x$size,
                       x$units, x$xscale, x$yscale)
    # Create id to distinguish separate symbols
    id <- rep(1:locn$n, each=length(x$object@x))
    # Generate grob representing symbols
    if (x$use.gc) {
        child <- do.call("polygonGrob",
                         c(list(x=locn$x, y=locn$y, id=id, 
                                default.units="inches",
                                gp=gpar(col=NA, fill=x$object@rgb,
                                        lineend=x$object@lineend,
                                        linejoin=x$object@linejoin,
                                        linemitre=x$object@linemitre)),
                           x$poly.args))
    } else {
        child <- do.call("polygonGrob",
                         c(list(x=locn$x, y=locn$y, id=id,
                                default.units="inches", vp=x$vp),
                           x$poly.args))
    }
    setChildren(x, gList(child))
}
    
setMethod("symbolize", signature(object="PictureStroke"),
          function(object,
                   x=unit(0.5, "npc"),
                   y=unit(0.5, "npc"),
                   size=unit(1, "npc"),
                   units="npc",
                   xscale=NULL, yscale=NULL, ..., use.gc=TRUE) {
              if (length(object@x) > 1) {
                  gTree(object=object, x=x, y=y, size=size,
                        units=units, xscale=xscale, yscale=yscale,
                        use.gc=use.gc, poly.args=list(...),
                        cl="symbolStroke")
              } else {
                  NULL
              }
          })

setMethod("symbolize", signature(object="PictureFill"),
          function(object,
                   x=unit(0.5, "npc"),
                   y=unit(0.5, "npc"),
                   size=unit(1, "npc"),
                   units="npc",
                   xscale=NULL, yscale=NULL, ..., use.gc=TRUE) {
              if (length(object@x) > 1) {
                  gTree(object=object, x=x, y=y, size=size,
                        units=units, xscale=xscale, yscale=yscale,
                        use.gc=use.gc, poly.args=list(...),
                        cl="symbolFill")
              } else {
                  NULL
              }
          })

setMethod("symbolize", signature(object="PictureText"),
          function(object,
                   x=unit(0.5, "npc"),
                   y=unit(0.5, "npc"),
                   size=unit(1, "npc"),
                   units="npc",
                   xscale=NULL, yscale=NULL, ..., use.gc=TRUE) {
              # Do not allow text in imported picture used for symbol (yet)
              NULL
          })

setMethod("symbolize", signature(object="PictureChar"),
          function(object,
                   x=unit(0.5, "npc"),
                   y=unit(0.5, "npc"),
                   size=unit(1, "npc"),
                   units="npc",
                   xscale=NULL, yscale=NULL, ..., use.gc=TRUE) {
              # Do not allow text in imported picture used for symbol (yet)
              NULL
          })

setMethod("symbolize", signature(object="Picture"),
          function(object, 
                   x=unit(0.5, "npc"),
                   y=unit(0.5, "npc"),
                   size=unit(1, "npc"),
                   units="npc",
                   ..., name=name, gp=gpar()) {
              gTree(children=do.call("gList",
                      lapply(object@paths, symbolize,
                             x, y, size, units,
                             xscale=object@summary@xscale,
                             yscale=object@summary@yscale,
                             ...)),
                    name=name, gp=gp)
          })

##################
# Draw entire picture
# or individual paths
##################
pictureGrob <- function(picture, 
                        x=0.5, y=0.5, width=1, height=1,
			just="centre", 
			exp=0.05, xscale=NULL, yscale=NULL,
                        distort=FALSE,
                        FUN=grobify, ...,
                        name=NULL, gp=gpar()) {
    grobify(picture, 
            x=x, y=y, width=width, height=height, just=just,
            exp=exp, xscale=xscale, yscale=yscale, distort=distort, FUN=FUN,
            ..., name=name, gp=gp)
}

grid.picture <- function(...) {
    grid.draw(pictureGrob(...))
}

symbolsGrob <- function(picture,
                        x=unit(0.5, "npc"),
                        y=unit(0.5, "npc"),
                        size=unit(1, "npc"),
                        units="npc",
                        ...,
                        name=NULL, gp=gpar()) {
    symbolize(picture,
              x=x, y=y, size=size, units=units,
              ..., name=name, gp=gp)
}

grid.symbols <- function(...) {
    grid.draw(symbolsGrob(...))
}

explodePaths <- function(picture) {
    picture@paths <- unlist(lapply(picture@paths, explode), recursive=FALSE)
    picture@summary <- new("PictureSummary",
                           numPaths=length(picture@paths),
                           xscale=picture@summary@xscale,
                           yscale=picture@summary@yscale)
    picture
}

picturePaths <- function(picture,
                         nr, nc,
                         col="black",
                         fill="light grey",
                         freeScales=FALSE, 
                         xscale=NULL, yscale=NULL,	      
                         label=function(n) { 
                             tg <- textGrob(n, x=0, y=0, 
                                            just=c("left", "bottom"),
                                            gp=gpar(fontsize=6))
                             grid.rect(x=0, y=0, height=unit(6, "points"),
                                       width=grobWidth(tg),
                                       just=c("left", "bottom"),
                                       gp=gpar(fill="white"))
                             grid.draw(tg) },
                         use.gc=TRUE) {
    if (missing(nr) || missing(nc)) {
      nrnc <- n2mfrow(picture@summary@numPaths)
      nr <- nrnc[1]
      nc <- nrnc[2]
    }
    if (is.null(xscale) || is.null(yscale)) {
        xscale <- picture@summary@xscale
	yscale <- picture@summary@yscale
    }
    if (freeScales) {
        pushViewport(viewport(layout=grid.layout(nr, nc)))        
    } else {
        pushViewport(viewport(layout=grid.layout(nr, nc,
                                widths=rep(diff(range(xscale)),
                                  nc),
                                heights=rep(diff(range(yscale)),
                                  nr),
                                respect=TRUE)))
    }
    for (i in 1:nr) {
        for (j in 1:nc) {
            pnum <- (i - 1)*nc + j
            if (pnum <= picture@summary@numPaths) {
                pushViewport(viewport(layout.pos.col=j,
                                      layout.pos.row=i))
                if (freeScales) {
                    pushViewport(viewport(width=0.95, height=0.95))
                    grid.rect(gp=gpar(col=col, fill=fill))
                    grid.picture(picture[pnum])
                    popViewport()
                } else {
                    pushViewport(viewport(width=0.95, height=0.95,
                                      xscale=xscale,
                                      yscale=yscale))
                    grid.rect(gp=gpar(col=col, fill=fill))
                    if (is.function(label)) {
                        label(pnum)
                    }
                    grob <- grobify(picture@paths[[pnum]], use.gc=use.gc)
                    if (is.null(grob)) {
                        grid.text("EMPTY\nPATH",
                                  gp=gpar(fontsize=6))
                    } else {
                        grid.draw(grob)
                    }
                    popViewport()
                }
                popViewport()
            }
        }
    }
    popViewport()
}

##################
# Retain old grobify() behaviour for back-compatibility

setGeneric("oldGrobify",
           function(object, ...) {
               standardGeneric("oldGrobify")
           })

picLinesGrob <- function(...) {
    grob(..., cl="picline")
}

drawDetails.picline <- function(x, recording) {
    # Figure out what lwd and lty really are
    lwd <- convertWidth(unit(x$lwd, "native"), "bigpts", valueOnly=TRUE)
    lty <- fixLTY(x$lty, x$lwd)
    grid.lines(x$x, x$y, x$default.units,
               gp=gpar(lwd=lwd, lty=lty, col=x$col))
}

# Individual path converted into grob
setMethod("oldGrobify", signature(object="PictureStroke"),
          function(object, ..., fillText, bgText, sizeByWidth, use.gc=TRUE) {
              if (length(object@x) > 1) {
                  if (use.gc) {
                      picLinesGrob(x=object@x, y=object@y,
                                   default.units="native",
                                   lwd=object@lwd,
                                   lty=object@lty,
                                   col=object@rgb, ...)
                  } else {
                      linesGrob(object@x, object@y,
                                default.units="native", ...)
                  }
              } else {
                  NULL
              }
          })

setMethod("oldGrobify", signature(object="PictureFill"),
          function(object, ..., fillText, bgText, sizeByWidth, use.gc=TRUE) {
              if (length(object@x) > 1) {
                  if (use.gc) {
                      polygonGrob(object@x, object@y, default.units="native",
                                  gp=gpar(col=NA, fill=object@rgb), ...)
                  } else {
                      polygonGrob(object@x, object@y,
                                  default.units="native", ...)
                  }
              } else {
                  NULL
              }
          })

setMethod("oldGrobify", signature(object="PictureText"),
          function(object, FUN=oldGrobify, ..., 
                   fillText=FALSE, bgText=.bgText.default,
                   sizeByWidth=TRUE, use.gc=TRUE) {
              if (length(object@letters) == 0) {
                  if (use.gc) {
                      pictureTextGrob(object@string,
                                      object@x, object@y,
                                      object@w, object@h,
                                      object@angle,
                                      object@letters,
                                      ...,
                                      sizeByWidth=sizeByWidth,
                                      gp=gpar(col=object@rgb))
                  } else {
                      pictureTextGrob(object@string,
                                      object@x, object@y,
                                      object@w, object@h,
                                      object@angle,
                                      object@letters,
                                      ...,
                                      sizeByWidth=sizeByWidth)
                  }
              } else {
                  gTree(string=as.character(object@string), 
                        children=do.call("gList",
                          lapply(object@letters, FUN=FUN, ...,
                                 fillText=fillText, bgText=bgText,
                                 sizeByWidth=sizeByWidth, use.gc=use.gc)),
                        cl="pictureletters")
              }
          })

setMethod("oldGrobify", signature(object="PictureChar"),
          function(object, FUN=oldGrobify, ..., 
                   fillText=FALSE, bgText=.bgText.default,
                   sizeByWidth=TRUE, use.gc=TRUE) {
              paths <- explode(object, fillText, bgText)
              do.call("gList", lapply(paths, FUN=FUN, ..., use.gc=use.gc))
          })
