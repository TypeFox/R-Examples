################################################################################
##
## $Id: exposure.R 346 2006-10-01 05:08:55Z enos $
##
## Class wrapping exposure data.
##
################################################################################

setMethod("show",
          signature(object = "exposure"),
          function(object){
            for(v in names(object@data)){
              cat(v, "\n")
              show(object@data[[v]])
              cat("\n")
            }
          }
          )


setMethod("summary",
          signature(object = "exposure"),
          function(object){
            show(object)
          }
          )

setMethod("plot",
          signature(x = "exposure", y = "missing"),
          function(x, ...){

            tlist <- list()
            for(var in names(x@data)){
              title <- paste(var, "exposure")
              data <- x@data[[var]]
              data$variable <- factor(data$variable, levels = rev(unique(data$variable)))
              tlist[[var]] <- barchart(variable ~ exposure, data = data, origin = 0, main = title)
            }

            .trellis.multiplot(tlist)
          }
          
          )
