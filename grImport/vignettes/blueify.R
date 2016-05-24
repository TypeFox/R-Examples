setGeneric("blueify",
           function(object, ...) {
               standardGeneric("blueify")
           })
setMethod("blueify", signature(object="PictureFill"),
          function (object, ...) {
              pathGrob(object@x, object@y, 
                       default.units="native",
                       gp=gpar(col=NA, 
                         fill=blueShade(object@rgb)), 
                       ...)
          })
setMethod("blueify", signature(object="PictureStroke"),
          function (object, ...) {
              polylineGrob(object@x, object@y, 
                           default.units="native",
                           gp=gpar(col=blueShade(object@rgb)), 
                           ...)
          })
