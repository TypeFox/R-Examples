setMethod("rowMeans",
          signature(x = "magpie"),
          function (x, na.rm = FALSE, dims = 1, ...) 
          {
            x <- rowMeans(as.array(x), na.rm=na.rm, dims=dims, ...)
            return(as.magpie(as.array(x)))
          }
          )