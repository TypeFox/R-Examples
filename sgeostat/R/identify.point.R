# A function to "identify" points in a "point" object...
# This function allows the user to identify points graphically.
assign("identify.point",
function(x,v="",...) {
    v <- x[[match(v,names(x))]]
    if(!is.null(v))
      identify(x$x,x$y,v,...)
    else
      identify(x$x,x$y,...)
})
