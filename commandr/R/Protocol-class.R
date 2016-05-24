## a protocol is a single step in a pipeline

## A protocol performs a stage in a particular way.
setClass("Protocol", contains = c("Command", "VIRTUAL"))

## get an instance of the stage to which this protocol belongs
setGeneric("stage", function(object, ...) standardGeneric("stage"))
## for the base classes
setMethod("stage", "Protocol",
          function(object, where = .externalCallerEnv()) {
            force(where)
            StageForProtocol(class(object), where)
          })

StageForProtocol <- function(name, where = topenv(parent.frame())) {
  if (!extends(name, "Protocol"))
    stop("Class '", name, "' is not a protocol class")
  ancestors <- names(getClass(name)@contains)
  protos <- sapply(c(name, ancestors), dequalifyProtocolName)
  stages <- sapply(findSubclasses("Stage", where), dequalifyStageName)
  stages <- stages[stages %in% protos]
  if (length(stages) == 1)
    Stage(stages[[1]])
  else if (length(stages) > 1)
    stop("Protocol '", name, "' inherits from multiple stages: ",
         paste("'", stages, "'", sep="", collapse=", "))
  else NULL
}



## name of specific method performed by a protocol
setGeneric("method", function(object, ...) standardGeneric("method"))

setMethod("method", "Protocol",
          function(object, where = .externalCallerEnv()) {
            force(where)
            role <- role(stage(object, where))
            decapitalize(sub(qualifyProtocolName(role), "", class(object)))
          })

setMethod("parameters", "Protocol", function(object) {
  ## simply return slots as a list
  slots <- lapply(slotNames(object),function(slot_name) slot(object, slot_name))
  names(slots) <- slotNames(object)
  slots
})

setMethod("show", "Protocol", function(object)
          {
            stage <- stage(object, .GlobalEnv)
            cat(displayName(stage), " (", displayName(object), ")\n", sep = "")
          })

### FIXME: should happen in setProtocol()
setMethod("pipeline", "Protocol", function(object)
          {
            if ("pipeline" %in% slotNames(object))
              object@pipeline
            else NULL
          })

## Get a protocol instance for the given stage and method
## eg Protocol("findPeaks", "matchedFilter") yields an instance of
## "ProtoFindPeaksMatchedFilter"

Protocol <- function(role, method = defaultMethod(role), ...)
{
  new(protocolClass(role, method), ...)
}

protocolClass <- function(role, method = NULL)
{
  if (is(role, "Stage"))
    role <- role(role)
  name <- paste(decapitalize(role), capitalize(method), sep="")
  qualifyProtocolName(name)
}

qualifyProtocolName <- function(name)
{
  paste("Proto", capitalize(name), sep = "")
}
dequalifyProtocolName <- function(name)
{
  decapitalize(sub("^Proto", "", name))
}

## Registration of protocols
setProtocol <- function(method, dispname = method, representation = list(),
                        fun, parent, prototype = list(), validity = NULL,
                        where = topenv(parent.frame()))
{
  ## no function, or "VIRTUAL" in 'parent', protocol is abstract
  virtual <- missing(fun) || "VIRTUAL" %in% parent
  parent <- setdiff(parent, "VIRTUAL")
  method <- decapitalize(method)
  ## resolve ancestors and find stage
  if (!extends(parent, "Protocol")){
    parent <- qualifyProtocolName(parent)
  }
  stage <- StageForProtocol(parent, where)
  if (is.null(stage))
    stop("Failed to derive a stage from parent class: '", parent, '"')
  stagename <- role(stage)
  ## class name directly computed from 'stage' and 'method'
  class <- protocolClass(stagename, method)
  if (dequalifyProtocolName(class) == stagename)
    stop("Protocol name conflicts with existing stage name '", stagename, "'")
  contains <- parent
  if (virtual)
    contains <- c(contains, "VIRTUAL")
  
  ## Transform representation to allow language objects (delayed evaluation)
  representation <- lapply(representation, function(cl) {
    union <- paste(cl, "language", sep="OR")
    if (!isClassUnion(union))
      setClassUnion(union, c(cl, "language"), where)
    union
  })


  ## add function formals to prototype
  if (!missing(fun)) {
    slots <- c(slotNames(parent), names(representation))
    params <- names(formals(fun)) %in% slots
    nonmissing <- params & nchar(sapply(formals(fun), deparse)) > 0
    prototype[names(formals(fun))[nonmissing]] <- formals(fun)[nonmissing]
  }
  ## create prototype without forcing argument evaluation
  prototype <- do.call("prototype", prototype, TRUE)
  setClass(class, representation, prototype, contains, validity, where = where)
  if (!missing(fun))
    setMethod("displayName", class, function(object) dispname, where = where)
  ## remember the 'stage' of the protocol
  ## creates a new instance since stages can be redefined
  ##setMethod("stage", class, function(object) Stage(stagename), where=where)
  ## remember the method name
  ##setMethod("method", class, function(object) method, where = where)
  if (!missing(fun)) {
    intype <- inType(stage)
    .fun <- fun
    if ("..." %in% names(formals(fun)))
      formal <- rep(TRUE, length(slotNames(class)))
    else formal <- slotNames(class) %in% names(formals(fun))
    setMethod("perform", c(class, intype),
              function(object, data, ...)
              {
                .data <- data
                slots <- parameters(object)                  
                expr <- quote({
                  nms <- sapply(names(formals(sys.function())), as.name)
                  do.call(.fun, c(list(.data), nms, list(...)))
                })
                result <- as.function(c(slots[formal], expr))()
                if (!is.null(result)) {
                  pipeline <- attr(data, "pipeline")
                  newPipe <- Pipeline(object)
                  if (is.null(pipeline))
                    pipeline <- newPipe
                  else pipeline <- c(pipeline, newPipe)
                  attr(result, "pipeline") <- pipeline
                  ##names(result@pipeline)[length(names(result@pipeline))] <- name(object)
                }
                result
              }, where = where)
    ## set a method with the same formals
    generic <- paste(stagename, method, sep=".")
    args <- formals(fun)
    parslots <- slots[!(slots %in% names(args))]
    args[parslots] <- attributes(getClass(parent)@prototype)[parslots]
    names(args)[1] <- "object"
    existing <- getGeneric(generic, where = where)
    def <- .dyngeneric(generic)
    if (is.null(existing) || !identical(formals(existing), formals(def)))
      setGeneric(generic, def, where = where)
    .proto <- list(stagename, method)
    .slots <- slots
    expr <- quote({
      mc <- as.list(match.call())[-c(1,2)]
      slots <- names(mc) %in% .slots
      protocol <- do.call("Protocol", c(.proto, mc[slots]))
      do.call("perform", c(list(protocol, object), mc[!slots]))
    })
    setMethod(generic, intype, as.function(c(args, expr)), where = where)
    if (is.null(defaultMethod(stage)))
      defaultMethod(stage) <- method
  }
  class
}


