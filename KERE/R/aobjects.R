setClass("kernel",representation("function",kpar = "list"))
setClass("kernelMatrix",representation("matrix"),prototype = structure(.Data =
                                                                         matrix()))

setClassUnion("listI", c("list","numeric","vector","integer","matrix"))
setClassUnion(
  "output", c(
    "matrix","factor","vector","logical","numeric","list","integer","NULL"
  )
)
setClassUnion("input", c("matrix","list"))
setClassUnion("kfunction", c("function","character"))
setClassUnion("mpinput", c("matrix","data.frame","missing"))
setClassUnion("lpinput", c("list","missing"))
setClassUnion("kpinput", c("kernelMatrix","missing"))



if (!isGeneric("kpar")) {
  if (is.function("kpar"))
    fun <- kpar
  else
    fun <- function(object)
      standardGeneric("kpar")
  setGeneric("kpar", fun)
}
setGeneric("kpar<-", function(x, value)
  standardGeneric("kpar<-"))
