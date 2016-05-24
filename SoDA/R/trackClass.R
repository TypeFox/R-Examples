setClass("track", representation(x ="numeric", y="numeric"))
setClass("track3", representation(z="numeric"), contains = "track")
setMethod("plot",
    signature("track", "missing"),
    function (x, y, ...) 
    {
        plot(x@x, x@y, asp = 1, ...)
    }
)
setMethod("plot", c("track3", "missing"), #version 2
          function(x, y, points = c(".", "o", "*"), ...) {
              which = as.integer(cut(x@z, length(points)))
              callNextMethod(x, pch = points[which], ...)
              })

setMethod("show", "track",
          function(object) {
              message("Object of class \"", class(object), "\"")
              print(cbind(x = object@x, y=object@y))
          })
