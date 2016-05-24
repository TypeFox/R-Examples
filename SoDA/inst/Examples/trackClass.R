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

## some example data
xy <- scan("Examples/gps1XY.txt", list(y = numeric(), x=numeric()))
tr = new("track", x=xy$x, y = xy$y)
t3 = new("track3", x = xy$x, y = xy$y, z=1:length(xy$x))
t3@x[1:10] = jitter(t3@x[1:10]) # to make comparisons interesting
