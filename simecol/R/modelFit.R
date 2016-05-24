############################################################
#  S4 class and methods for model fits                     #
############################################################


setClass("modelFit",
         representation(
           par = "numeric",
           message = "ANY",
           value = "ANY",
           list = "ANY"
         ))

setMethod("coef", "modelFit",
          function(object, ...) {
            object@par
          }
)

setMethod("deviance", "modelFit",
          function(object, ...) {
            object@value
          }
) 

setMethod("summary", "modelFit",
          function(object, ...) {
            cat("Model parameters", object@par, "\n")
            cat("Deviance:", object@deviance, "\n")
            cat("Message:", object@message, "\n")
          })


## accessor methods for S3 compatibility (read-only)
setMethod("$", "modelFit",
          function(x, name)
          {
            ## 'name' is a character(1)
            slot(x, name)
          })


setMethod("[", "modelFit",
          function(x, i, j, ..., drop=TRUE)
          {
            x@list[i, j, ..., drop=drop]
          })


setMethod("[[", "modelFit",
          function(x, i, j, ...)
          {
            x@list[[i, j, ...]]
          })

