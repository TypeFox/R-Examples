setMethod("show",
          signature(object = "SpatialStreamNetwork"),
          function (object) {
              cat("Object of class Spatial Stream Network\n\n")

              nobs <- dim(object@obspoints@SSNPoints[[1]]@point.data)
              nobs <- matrix(nobs,1,)
              np <- length(object@predpoints@SSNPoints)
              if(np > 0) {
                  for(i in 1:np) {
                      nobs <- rbind(nobs,dim(object@predpoints@SSNPoints[[i]]@point.data))
                  }
              }

              cat("Object includes observations on",nobs[1,2],"variables across",nobs[1,1],"sites within the bounding box\n")
              print(bbox(object))
              cat("\n")

              if(nrow(nobs)==2) {
                  cat("Object also includes", nrow(nobs)-1, "set of prediction points with",sum(nobs[,1]) - nobs[1,1],"locations\n\n")
              } else if(nrow(nobs)>2) {
                  cat("Object also includes", nrow(nobs)-1, "sets of prediction points with a total of",sum(nobs[,1]) - nobs[1,1],"locations\n\n")
              }

              cat("Variables recorded are (found using names(object)):\n")
              print(names(object))

              cat("Generic functions that work with this object include names, plot, print, summary, hist, boxplot and qqnorm\n")
          })
