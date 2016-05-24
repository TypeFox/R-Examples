setGeneric("grayify",
           function(object, ..., min=0) {
               standardGeneric("grayify")
           })
setMethod("grayify", signature(object="PictureStroke"),
          function(object, ..., min=0) {
              linesGrob(object@x, object@y, default.units="native",
                        gp=gpar(col=gray(min)), ...)
          })
setMethod("grayify", signature(object="PictureFill"),
          function(object, ..., min=0) {
              rgbcol <- hex2RGB(object@rgb)
              luvcol <- as(rgbcol, "LUV")
              # Extract just Luminance component as gray level
              # Rescale to [0.75, 1] to lighten "black" to "gray75"
              graycol <- min + coords(luvcol)[,1]*(1 - min)/100
              polygonGrob(object@x, object@y, default.units="native",
                          gp=gpar(lwd=object@lwd,
                            col=NA, fill=gray(graycol)), ...)
          })
