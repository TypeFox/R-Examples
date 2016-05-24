setClass("irregTS", representation(time = "DateTime"), 
   contains = "structure")

setMethod("[", "irregTS", function (x, i, j, ..., drop = TRUE)
          new("irregTS", x@.Data[i], time = x@time[i])
          )

setMethod("show", "irregTS", function(object) {
         cat("Object of class \"", class(object), "\"\n")
          tt <- object@time
          start = min(tt, na.rm = TRUE)
          cat("Times starting at ", format(start), "\n")
         x = object@.Data
         names(x) = tt - start
         show(x)
     }
          )
