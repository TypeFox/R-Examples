# Conversions between class

setAs("loca.p", "data.frame",
          function(from, to)
             data.frame("w" = from@w, "x" = from@x, "y" = from@y)
      )

setAs("data.frame", "loca.p", 
          function(from, to)
            new("loca.p", w = from$w, x = from$x, y = from$y)
      )


setAs("loca.p", "matrix",
          function(from, to)
            cbind(from@x, from@y, from@w)       
      )

setAs("matrix", "loca.p",
          function(from, to)
            {
              if (dim(from)[2] == 2) loca.p(x=from[,1], y=from[,2])
              else if (dim(from)[2] == 3) loca.p(x=from[,1], y=from[,2], w=from[,3])
              else stop(gettext("The second dimension of matrix must be 2 or 3.", domain = "R-orloca"))
            }
      )

#
# The following is por S3 compatibility, mainly for documentation check
#
setGeneric("as.loca.p", function(x, ...) standardGeneric("as.loca.p"))
as.loca.p.matrix <- function(x,...) as(x, "loca.p")
setMethod("as.loca.p", "matrix", as.loca.p.matrix)
as.loca.p.data.frame <- function(x, ...) as(x, "loca.p")
setMethod("as.loca.p", "data.frame", as.loca.p.data.frame)
as.matrix.loca.p <- function(x, ...) as(x, "matrix")
setMethod("as.matrix", "loca.p", as.matrix.loca.p)
as.data.frame.loca.p <- function(x, row.names = NULL, optional = FALSE, ...) as(x, "data.frame")
setMethod("as.data.frame", "loca.p", as.data.frame.loca.p)
