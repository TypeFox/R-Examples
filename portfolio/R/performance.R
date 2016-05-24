################################################################################
##
## $Id: performance.R 393 2007-01-23 15:10:39Z enos $
##
## Class wrapping performance data.
##
################################################################################

setMethod("initialize",
          signature(.Object = "performance"),
          function(.Object, ...){
            if(nrow(.Object@ret.detail) == 0){
              row.names(.Object@ret.detail) <- integer(0)
            }
            .Object
          }
          )

setMethod("show",
          signature(object = "performance"),
          function(object){
            if(length(object@ret) > 0){

              ret <- object@ret
              ret.tag  <- ifelse(abs(ret) > 0.01, "%", "bps")
              ret <- ifelse(ret.tag == "%", ret * 100, ret * 100 * 100)
              ret <- round(ret, digits = 2)
              
              cat(paste("Total return: ", ret, ret.tag, "\n\n"))
              if(nrow(object@ret.detail) > 0){
                x <- object@ret.detail
                x <- x[order(x$contrib, na.last = NA),]
                cat("Best/Worst performers:\n")

                if(nrow(x) < 10){
                  show(x)
                }
                else{
                  show(rbind(head(x, n = 5),
                             tail(x, n = 5)))
                }
              }
            }
            else{
              cat("Object of class performance with no return data.\n")
            }
          }
          )

setMethod("summary",
          signature(object = "performance"),
          function(object){
            show(object)
          }
          )

setMethod("plot",
          signature(x = "performance", y = "missing"),
          function(x){
            ret <- x@ret
            ret.tag  <- ifelse(abs(ret) > 0.01, "%", "bps")
            ret <- ifelse(ret.tag == "%", ret * 100, ret * 100 * 100)
            ret <- round(ret, digits = 2)

            if(nrow(x@ret.detail) > 0){
              y <- x@ret.detail
              y <- y[order(y$contrib, na.last = NA),]
              y$id <- factor(y$id, levels = unique(y$id))
              if(nrow(y) > 10){
                y <- rbind(head(y, n = 5),
                           tail(y, n = 5))
              }

              print(barchart(id ~ contrib, data = y, origin = 0,
                             main = paste("Performance: ", ret, ret.tag, sep = "")))
                             
            }
            else{
              stop("Nothing to plot")
            }
          }
          )

