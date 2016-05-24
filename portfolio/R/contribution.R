################################################################################
##
## $Id: contribution.R 346 2006-10-01 05:08:55Z enos $
##
## Class wrapping contribution data.
##
################################################################################

setMethod("show",
          signature(object = "contribution"),
          function(object){
            for(v in names(object@data)){
              cat(v, "\n")
              show(object@data[[v]])
              cat("\n")
            }
          }
          )

setMethod("summary",
          signature(object = "contribution"),
          function(object){
            show(object)
          }
          )

setMethod("plot",
          signature(x = "contribution", y = "missing"),
          function(x, ...){

            tlist <- list()
            for(var in names(x@data)){
              title <- paste("roic by", var)
              data <- x@data[[var]]
              data$variable <- factor(data$variable, levels = rev(unique(data$variable)))
              tlist[[var]] <- barchart(variable ~ roic, data = data, origin = 0, main = title)
            }
            
            .trellis.multiplot(tlist)
          }
          )
