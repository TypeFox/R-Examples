#########################################################################/**
# @set class=default
# @RdocMethod getAttributes
# @aliasmethod getAttribute
# @aliasmethod hasAttribute
# @aliasmethod setAttributes
# @aliasmethod setAttribute
# @aliasmethod getMetadata
# @aliasmethod setMetadata
#
# @title "Gets and sets attributes of an object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{object}{An object.}
#   \item{private}{If @TRUE, attributes starting with a period are
#         also returned, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list, @NULL or a modified object itself.
# }
#
# @author
#
# @keyword internal
#*/#########################################################################
setMethodS3("getAttributes", "default", function(object, private=FALSE, ...) {
  attrs <- attributes(object)
  keys <- names(attrs)
  keys <- setdiff(keys, c("class", "names"))

  # Exclude private attributes?
  if (!private) {
    pattern <- sprintf("^[%s]", paste(c(base::letters, base::LETTERS), collapse=""))
    keys <- keys[regexpr(pattern, keys) != -1L]
  }

  attrs <- attrs[keys]
  attrs
})

setMethodS3("getAttribute", "default", function(object, name, default=NULL, private=TRUE, ...) {
  attrs <- getAttributes(object, private=private, ...)
  if (!is.element(name, names(attrs))) {
    attr <- default
  } else {
    attr <- attrs[[name]]
  }
  attr
})

setMethodS3("hasAttribute", "default", function(object, name, private=TRUE, ...) {
  attrs <- getAttributes(object, private=private, ...)
  is.element(name, names(attrs))
})

setMethodS3("setAttributes", "default", function(object, attrs, ...) {
  # Argument 'attrs':
  if (is.null(attrs)) {
    return(invisible(object))
  }
  if (!is.list(attrs)) {
    throw("Cannot set attributes. Argument 'attrs' is not a list: ", mode(attrs)[1L])
  }


  # Current attributes
  attrsD <- attributes(object)

  # Update/add new attributes
  keys <- names(attrs)
  keys <- setdiff(keys, c("class", "names"))
  for (key in keys) {
    attrsD[[key]] <- attrs[[key]]
  }

  attributes(object) <- attrsD

  invisible(object)
})

setMethodS3("setAttribute", "default", function(object, name, value, ...) {
  attrs <- list(value)
  names(attrs) <- name
  setAttributes(object, attrs, ...)
})




setMethodS3("getMetadata", "default", function(object, name=NULL, default=NULL, local=FALSE, ...) {
  res <- getAttribute(object, "metadata", default=list())
  if (!local) {
    isLocal <- is.element(names(res), "source")
    res <- res[!isLocal]
  }
  if (!is.null(name)) {
    if (is.element(name, names(res))) {
      res <- res[[name]]
    } else {
      res <- default
    }
  }
  res
}, protected=TRUE)


setMethodS3("setMetadata", "default", function(object, metadata=NULL, name, value, ...) {
  data <- getMetadata(object, local=TRUE)

  if (!is.null(metadata)) {
    for (name in names(metadata)) {
      data[[name]] <- metadata[[name]]
    }
  } else {
    data[[name]] <- value
  }

  setAttribute(object, "metadata", data)
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2015-02-04
# o Now all attribute and metadata functions are for "default" objects.
# o Created from RspNnn.R files.
##############################################################################
