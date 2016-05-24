
# Draw a "picture" object in the current plot region

setGeneric("drawPath",
           function(p, trans, ...) {
               standardGeneric("drawPath")
           })

setMethod("drawPath", signature(p="PictureStroke"),
          function(p, trans, ...) {
              lwd <- (trans(p@lwd, 0)$x - trans(0, 0)$x)/xinch()*72
              lty <- fixLTY(p@lty, p@lwd)
              paths <- explode(p)
              if (is.list(paths)) {
                  pathX <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@x) }))[-1]
                  pathY <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@y) }))[-1]
              } else {
                  pathX <- paths@x
                  pathY <- paths@y
              }
              lines(trans(pathX, pathY), col=p@rgb, lwd=lwd, lty=lty)
          })

setMethod("drawPath", signature(p="PictureFill"),
          function(p, trans, ...) {
              paths <- explode(p)
              if (is.list(paths)) {
                  pathX <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@x) }))[-1]
                  pathY <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@y) }))[-1]
              } else {
                  pathX <- paths@x
                  pathY <- paths@y
              }
              polypath(trans(pathX, pathY),
                       rule=switch(p@rule,
                         nonzero="winding", "evenodd"),
                       col=p@rgb, border=NA)
          })

setMethod("drawPath", signature(p="PictureText"),
          function(p, trans, ...) {
              if (length(p@letters) == 0) {
                  text(trans(p@x, p@y), labels=p@string, col=p@rgb,
                       srt=p@angle)
              } else {
                  lapply(p@letters, drawPath, trans, ...)
              }
          })

setMethod("drawPath", signature(p="PictureChar"),
          function(p, trans, ...) {
              paths <- explode(p, FALSE, NA)
              if (is.list(paths)) {
                  pathX <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@x) }))[-1]
                  pathY <- unlist(lapply(paths,
                                         function(pp) { c(NA, pp@y) }))[-1]
              } else {
                  pathX <- paths@x
                  pathY <- paths@y
              }
              polypath(trans(pathX, pathY), 
                       rule="winding",
                       col=p@rgb, border=NA)
          })

picture <- function(picture, xleft, ybottom, xright, ytop, ...) {
    # Determine the conversion from picture coordinates
    # to plot coordinates
    dx <- xright - xleft
    dy <- ytop - ybottom
    xscale <- picture@summary@xscale
    yscale <- picture@summary@yscale
    trans <- function(x, y) {
        list(x=xleft + (x - min(xscale))/diff(xscale) * dx,
             y=ybottom + (y - min(yscale))/diff(yscale) * dy)
    }
    invisible(lapply(picture@paths, drawPath, trans))
}
