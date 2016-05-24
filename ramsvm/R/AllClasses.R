methods::setClass(Class = "ramsvm",
                  contains = "VIRTUAL")

methods::setClass("linear_RAM", 
                  slots = c("x" = "matrix",
                            "y" = "vector",
                            "y.name" = "vector",
                            "k" = "integer",
                            "gamma" = "numeric",
                            "weight" = "vector",
                            "lambda" = "vector",
                            "beta" = "list",
                            "beta0" = "list",
                            "epsilon" = "numeric",
                            "warm" = "matrix",
                            "call" = "call"),
                  contains = c("ramsvm") )

methods::setClass("kernel_RAM", 
                  slots = c("x" = "matrix",
                            "y" = "vector",
                            "y.name" = "vector",
                            "k" = "integer",
                            "gamma" = "numeric",
                            "weight" = "vector",
                            "lambda" = "vector",
                            "kernel" = "character",
                            "kparam" = "numeric",
                            "beta" = "list",
                            "beta0" = "list",
                            "epsilon" = "numeric",
                            "warm" = "matrix",
                            "call" = "call"),
                  contains = c("ramsvm") )


methods::setGeneric(name = "Pred", 
                    def = function(object, ...){
                            standardGeneric("Pred")
                          })

methods::setMethod(f = "Pred", 
                   signature = c(object = "linear_RAM"), 
                   definition = function(object, x, beta, beta0, ...) {

                                  res <- pred_linear(x.test = x,
                                                     beta = beta,
                                                     beta0 = beta0,
                                                     k = object@k)

                                  return( res )
                                } )

methods::setMethod(f = "Pred", 
                   signature = c(object = "kernel_RAM"), 
                   definition = function(object, x, beta, beta0, ...){

                                  res <- pred_kernel(x.test = x,
                                                     x = object@x,
                                                     kernel = object@kernel,
                                                     kparam = object@kparam,
                                                     beta = beta,
                                                     beta0 = beta0,
                                                     k = object@k)

                                  return( res )
                                } )

