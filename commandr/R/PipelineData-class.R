## PipelineData: We define an S4 class with a slot/attribute
## "pipeline" for storing the provenance of the object. However, these
## methods apply to any object with an attribute of the same name.

setClass("PipelineData", representation(pipeline = "Pipeline"))

setMethod("pipeline", "ANY",
          function(object, ancestry = TRUE, local = TRUE)
          {
            pipeline <- attr(object, "pipeline")
            locals <- pipeline@.Data
            me <- sapply(sapply(pipeline, outType), extends, class(object))
            if (any(!me))
              locals <- tail(pipeline, -tail(which(!me), 1))
            ancestors <- list()
            if (ancestry)
              ancestors <- head(pipeline, -length(locals))
            pipeline@.Data <- c(ancestors, if (local) locals)
            pipeline
          })

## explore the data in the context of the last applied protocol
setMethod("explore", c("ANY", "missing"),
          function(object, protocol, ...)
          {
            proto <- NULL
            pipeline <- pipeline(object)
            if (length(pipeline))
              proto <- tail(pipeline, 1)[[1]]
            explore(object, proto, ...)
          })
