### The pipeline is a series of protocols that implement stages
setClassUnion("OptionalCharacter", c("character", "NULL"))

setClass("Pipeline", representation(displayName = "OptionalCharacter"),
         contains = "list")

Pipeline <- function(..., displayName = NULL) {
  protos <- list(...)
  chars <- sapply(protos, is.character)
  protos[chars] <- lapply(protos[chars], Protocol)
  if (!all(sapply(protos, is, "Protocol")))
    stop("All arguments in '...' must be Protocol instances")
  if (!is.null(displayName) &&
      (!is.character(displayName) || length(displayName) != 1))
    stop("'displayName' should be a single character string or NULL")
  new("Pipeline", protos, displayName = displayName)
}

## protocol accessors
setGeneric("protocol", function(object, ...) standardGeneric("protocol"))
setGeneric("protocol<-", function(object, ..., value)
           standardGeneric("protocol<-"))

## extract a pipeline from an object
setGeneric("pipeline", function(object, ...) standardGeneric("pipeline"))

## The display name accessor
setMethod("displayName", "Pipeline", function(object) {
  object@displayName
})

## name of type to accept as input
## this is a single class, but it could be a class union
setGeneric("inType", function(object, ...) standardGeneric("inType"))

setMethod("inType", "Pipeline", function(object) {
  first <- head(object, 1)
  if (length(first))
    inType(first[[1]])
  else NULL
})

## name of type to produce
setGeneric("outType", function(object, ...) standardGeneric("outType"))

setMethod("outType", "Pipeline", function(object) {
  last <- tail(object, 1)
  if (length(last))
    outType(last[[1]])
  else NULL
})

## parameters controlling protocol behavior
setGeneric("parameters", function(object) standardGeneric("parameters"))

setMethod("parameters", "Pipeline", function(object) {
  lapply(object, parameters)
})

## return the longest range of protocols that goes from inType to outType
setGeneric("pipeline", function(object, ...) standardGeneric("pipeline"))
setMethod("pipeline", "Pipeline",
          function(object, intype = "ANY", outtype = "ANY")
          {
            inmatch <- sapply(sapply(object, inType), extends, intype)
            outmatch <- sapply(sapply(object, outType), extends, outtype)
            pipeline <- initialize(object, list())
            if (any(inmatch) && any(outmatch)) {
              pipeline <- object
              pipeline@.Data <-
                pipeline[which(inmatch)[1]:tail(which(outmatch),1)]
            }
            pipeline
          })

setGeneric("findProtocols", function(object, ...)
           standardGeneric("findProtocols"))
setMethod("findProtocols", "Pipeline",
          function(object, role, method = character())
          {
            which(sapply(object, is, protocolClass(role, method)))
          })

setMethod("protocol", "Pipeline",
          function(object, role, method = character())
          {
            protos <- findProtocols(object, role, method)
            if (!length(protos))
              NULL
            else object[[protos[1]]]
          })

setReplaceMethod("protocol", "Pipeline",
                 function(object, role, value)
                 {
                   protos <- findProtocols(object, role)
                   if (length(protos))
                     object@.Data[[protos[1]]] <- value
                   else object@.Data <- c(object, value)
                   object
                 })

## Get Pipelines up to a certain protocol or outtype
setMethod("head", "Pipeline",
          function(x, n = 6L, role, method = character(), outtype)
          {
            if (!missing(outtype)) {
              outmatch <- sapply(sapply(x, outType), extends, outtype)
              if (!any(outmatch))
                return(NULL)
              n <- which(outmatch)[1]
            } else if (!missing(role)) {
              protos <- findProtocols(x, role, method)
              if (!length(protos))
                return(NULL)
              else n <- protos[1]
            }
            initialize(x, callNextMethod(x, n = n))
          })

## Or starting at a protocol
setMethod("tail", "Pipeline",
          function(x, n = 6L, role, method = character(), intype)
          {
            if (!missing(intype)) {
              inmatch <- sapply(sapply(x, inType), extends, intype)
              if (!any(inmatch))
                return(NULL)
              n <- tail(which(inmatch), 1)
            } else if (!missing(role)) {
              protos <- findProtocols(x, role, method)
              if (!length(protos))
                return(NULL)
              else n <- tail(protos, 1) - length(x) - 1L
            }
            initialize(x, callNextMethod(x, n = n))
          })

## perform an operation on a data structure and return the result
## methods are defined via setProtocol()
setGeneric("perform",
           function(object, data = NULL, ...) standardGeneric("perform"))

## Perform all component protocols
setMethod("perform", "Pipeline", function(object, data)
          {
            for (proto in object)
              data <- perform(proto, data)
            data
          })

setMethod("show", "Pipeline", function(object) {
  cat("Pipeline with", length(object@.Data), "protocol(s):\n")
  for (i in seq_along(object@.Data)) {
    cat("[", i, "] ", sep = "")
    show(object[[i]])
  }
})

## NOTE: for all primitives, callNextMethod() will return a plain list
## need validation and consider inheritance
.valid.IOtype <- function(x){
  if(length(x) > 1){
  itype <- sapply(x, inType)
  otype <- sapply(x, outType)
  N <- length(itype)
  idx <- sapply(1:(N-1), function(i){
    extends(otype[i], itype[i+1])
  })
  id <- which(!idx)
  if(length(id)){
    msg <- character()
    for(i in id){
      msg <- c(msg, paste0("Pipeline ", i, " (outtype): ", otype[i], 
                          "' doesn't match Pipeline ", i +1, " (intype): ",
                           itype[i + 1]))
    }
    return(msg)
  }else{
    return(TRUE)
  }}else{
    return(TRUE)
  }
}


setValidity("Pipeline", .valid.IOtype)

setMethod("c", "Pipeline", function(x, ..., recursive = FALSE) {
  initialize(x, callNextMethod())
})

setMethod("[", "Pipeline",
          function(x, i, j, ..., drop=FALSE) {
            initialize(x, callNextMethod())
          })

setReplaceMethod("[", "Pipeline",
                 function(x, i, j, ..., value) {
                   ## FIXME: commented code doesn't work properly
                  ## initialize(x, callNextMethod())
                   x@.Data[i] <- value
                   validObject(x)
                   x
                 })
