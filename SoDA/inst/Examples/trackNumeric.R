setClass("trackNumeric",
         contains = c("numeric", "track"))

setValidity("trackNumeric",  function(object) {
             slotLengths <- c(length(object@.Data), length(object@x),
                              length(object@y))
             if(length(unique(slotLengths)) == 1)
               TRUE
             else
               paste("Differing lengths for data, x, y:",
                     paste(slotLengths, collapse = ", "))
         })

setMethod("show", "trackNumeric",
          function(object) {
              message("Object of class: \"", class(object), "\"")
              show(cbind(object@.Data, x = object@x, y = object@y))
          })

setMethod("plot", c("trackNumeric", "missing"), #version 2
          function(x, y, points = c(".", "o", "*"), ...) {
              which = as.integer(cut(x@.Data, length(points)))
              callNextMethod(x, pch = points[which], ...)
              })
setMethod("[",
    signature("trackNumeric", "numeric"),
    function (x, i, j, ..., drop = TRUE) 
    {
        if(all(i>0)) {
        	## must be monotone so it keeps the "track" concept
        	diffi <- diff(i)
        	if(all(diffi>0) || all(diffi<0)) {} #OK
        	else
        	  stop("Only subsets that extract in a single direction produce a valid track object")
        	}
        	x@.Data <- x@.Data[i]
        	x@x <- x@x[i]
        	x@y <- x@y[i]
        	x
    }
)
