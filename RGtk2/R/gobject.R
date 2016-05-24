##### GObject wrapping #####

# GType support

# currently just converts type name to GType object, if x isn't one already
as.GType <- function(x)
{
  mapping <- c("integer" = "gint", "character" = "gchararray", "logical" = "gboolean",
    "numeric" = "gdouble", "raw" = "guchar", "externalptr" = "gpointer") 
  type <- x
  if (is.character(type)) {
    if (type %in% names(mapping))
      type <- mapping[[type]]
    type <- try(gTypeFromName(type), TRUE)
    if (inherits(type, "try-error")) {
      func <- paste(tolower(substring(x, 1, 1)), substring(x, 2), "GetType", sep="")
      if (exists(func))
        type <- do.call(func, list())
      }
  }
  if (!inherits(type, "GType"))
    stop("Cannot convert ", x, " to GType")
  type
}

interface <-
function(obj)
{
 attr(obj, "interfaces")
}

gTypeGetAncestors <-
function(type)
{
  type <- as.GType(type)
  .Call("R_getGTypeAncestors", type, PACKAGE = "RGtk2")
}

gTypeGetInterfaces <-
function(type)
{
  type <- as.GType(type)
  .Call("R_getInterfaces", type, PACKAGE = "RGtk2")
}

gTypeGetClass <-
function(type)
{
  type <- as.GType(type)
  ancestors <- gTypeGetAncestors(type)
  class_ptr <- .Call("R_getGTypeClass", type, PACKAGE = "RGtk2")
  class(class_ptr) <- c(paste(ancestors, "Class", sep=""), class(class_ptr))
  class_ptr
}

gTypeFromName <-
function(name)
{
 .Call("R_gTypeFromName", as.character(name), PACKAGE = "RGtk2")
}

print.GType <- function(x, ...) {
  cat("GType identifier for '", attr(x, "name"), "'\n", sep = "")
}

# GSignal support

GSignalFlags <- c(
  "run-first"  = 1,
  "run-last"  = 2,
  "run-cleanup"  = 4,
  "no-recurse"  = 8,
  "detailed"  = 16,
  "action"  = 32,
  "no-hooks"  = 64
)

GConnectFlags <- c(
  "after" = 1,
  "swapped" = 2
)

connectSignal <- gSignalConnect <-
function(obj, signal, f, data = NULL, after = FALSE, user.data.first = FALSE)
{
  useData <- missing(data) == FALSE
  checkPtrType(obj, "GObject")

  if(is.null(f))
    stop("You've specified NULL as the action in setting a callback. Did you mean to use quote()")

  if(is.expression(f)) {
    f <- f[[1]]
  }

  if(!( is.expression(f) || is.function(f) || is.call(f))) {
    stop(paste("Callback action must be an expression, a call or a function, but instead is of type", typeof(f), ". Did you forget to use quote()"))
  }

  invisible(.Call("R_connectGSignalHandler", obj, f, as.character(signal), data,
                  useData, as.logical(after), as.logical(user.data.first),
                  PACKAGE = "RGtk2"))
}

print.CallbackID <- function(x, ...)
  cat("Connection to '", names(x), "': ", x, "\n", sep = "")

gSignalHandlerDisconnect <-
function(obj, id)
{
 checkPtrType(obj, "GObject")
 .Call("R_disconnectGSignalHandler", obj, as.integer(id), PACKAGE = "RGtk2")
}

gSignalHandlerBlock <-
function(obj, id)
{
  checkPtrType(obj, "GObject")
 .Call("R_blockGSignalHandler", obj, as.integer(id), TRUE, PACKAGE = "RGtk2")
}

gSignalHandlerUnblock <-
function(obj, id)
{
  checkPtrType(obj, "GObject")
 .Call("R_blockGSignalHandler", obj, as.integer(id), FALSE, PACKAGE = "RGtk2")
}

gSignalStopEmission <-
function(obj, signal, detail = NULL)
{
  if (!is.null(detail))
    signal <- paste(signal, detail, sep="::")
  .Call("R_gSignalStopEmission", obj, signal, PACKAGE = "RGtk2")
}

gObjectGetSignals <-
function(obj)
{
  checkPtrType(obj, "GObject")
  type <- class(obj)[1]
  els <- gTypeGetSignals(type)
  els
}

gTypeGetSignals <-
function(type)
{
  if(is.character(type))
    type <- as.GType(type)

  checkPtrType(type, "GType")
  els <- .Call("R_getGSignalIdsByType", type, PACKAGE = "RGtk2")

  names(els) <- sapply(els, function(x) names(x))

  els
}


gSignalGetInfo <-
function(sig)
{
 checkPtrType(sig, "GSignalId")
 .Call("R_getGSignalInfo", sig, PACKAGE = "RGtk2")
}

gSignalEmit <-
function(obj, signal, ..., detail = NULL)
{
  checkPtrType(obj, "GObject")
  args <- list(...)
  signal <- as.character(signal)
  if (!is.null(detail))
    signal <- paste(signal, detail, sep="::")
  .RGtkCall("R_gSignalEmit", obj, signal, args, PACKAGE = "RGtk2")
}

# GObject properties

names.GObject <-
  #
  # return a vector of the names of the properties
  # available for the given GObject, collapsing over
  # all the inherited classes and removing the class::
  # prefix.
  #
function(x)
{
  names(gObjectGetPropInfo(x, parents = TRUE, collapse = TRUE))
}

gObjectGetPropInfo <-
function(obj, parents = TRUE, collapse = TRUE)
{
  checkPtrType(obj, "GObject")
  real_classes <- class(obj)[-length(class(obj))]
  props <- lapply(real_classes, gTypeGetPropInfo)
  if (parents && collapse)
    return(props[[1]])
  # props is a list containing the properties for each class in the hierarchy
  # as well as the parents of that class. We must remove the duplicates.
  n_dups <- c(sapply(props, length), 0)
  stripped <- lapply(1:length(props), function(ind) 
    if (n_dups[ind+1] > 0)
      props[[ind]][-(1:n_dups[ind+1])]
    else props[[ind]])
  names(stripped) <- real_classes
  result <- stripped
  if (!parents)
    result <- stripped[[1]]
  result
}

gTypeGetPropInfo <-
function(type)
{
  type <- as.GType(type)
  if (!("GObject" %in% gTypeGetAncestors(type)))
    stop("Cannot retrieve properties, because type is not a GObject type")
  
  .RGtkCall("R_getGTypeParamSpecs", type)
}

gObjectGet <-
function(obj, ..., drop = T)
{
   checkPtrType(obj, "GObject")
   props <- .Call("R_getGObjectProps", obj, as.character(c(...)), PACKAGE = "RGtk2")
   if (drop && length(props) == 1)
     props[[1]]
   else props
}

"[.GObject" <-
function(obj, value, ...)
{
 gObjectGet(obj, c(value, ...))
}

gObjectSet <-
function(obj, ...)
{
  args <- list(...)
  checkPtrType(obj, "GObject")
  if(any(names(args) == ""))
    stop("All values must have a name")

  invisible(.RGtkCall("R_setGObjectProps", obj, args, PACKAGE = "RGtk2"))
}

"[<-.GObject" <-
function(obj, propNames, value)
{
  value <- list(value)
  names(value) <- propNames
  .RGtkCall("R_setGObjectProps", obj, value, PACKAGE = "RGtk2")
  obj
}

gObject <- gObjectNew <-
function(type, ...)
{
  args <- list(...)
  type <- as.GType(type)
  if (!("GObject" %in% gTypeGetAncestors(type)))
    stop("GType must inherit from GObject")
  if(any(names(args) == ""))
    stop("All values must have a name")

  invisible(.RGtkCall("R_gObjectNew", type, args, PACKAGE = "RGtk2"))
}

## Parameter specifications

GParamFlags <- c("readable" = 1, "writable" = 2, "construct" = 4, 
  "construct-only" = 8, "lax-validation" = 16, "static-name" = 32,
  "private" = 32, "static-nick" = 64, "static-blurb" = 128)
  
gParamSpec <-
function(type, name, nick = NULL, blurb = NULL, flags = NULL, ...)
{
  # map type to param spec type, pass on the args
  
  spec <- list(name = name, nick = nick, blurb = blurb, flags = flags, ...)
  
  if (type == "integer")
    param_type <- "GParamInt"
  else if (type == "numeric")
    param_type <- "GParamDouble"
  else if (type == "logical")
    param_type <- "GParamBoolean"
  else if (type == "character")
    param_type <- "GParamString"
  else if (type == "raw")
    param_type <- "GParamUChar"
  else if (type == "R")
    param_type <- "RGtkParamSexp"
  else param_type <- type
  
  class(spec) <- c(param_type)
  
  as.GParamSpec(spec)
}

as.GParamSpec <- 
function(x)
{
  type <- sub(".*Param", "", class(x)[1])
  
  fields <- NULL
  common_fields <- c("name", "nick", "blurb", "flags")
  if (type %in% c("Boolean", "String", "Unichar"))
    fields <- "default.value"
  else if (type == "Flags")
    fields <- c("flags.type", "default.value")
  else if (type == "Enum")
    fields <- c("enum.type", "default.value")
  else if (type %in% c("Char", "UChar", "Int", "UInt", "ULong", "Long", "UInt64",
      "Int64", "Float", "Double"))
    fields <- c("minimum", "maximum", "default.value")
  else if (type == "Param")
    fields <- "param.type"
  else if (type == "Boxed")
    fields <- "boxed.type"
  else if (type == "Object")
    fields <- "object.type"
  else if (type == "ValueArray")
    fields <- "element.spec"
  else if (type == "GType")
    fields <- "is.a.type"
  else if (type == "Sexp")
    fields <- c("s.type", "default.value")
  
  x <- as.struct(x, c(class(x)[1], "GParamSpec"), c(common_fields, fields))

  x[[1]] <- as.character(x[[1]])
  x[[2]] <- as.character(x[[2]])
  x[[3]] <- as.character(x[[3]])
  
  if (is.null(x[[4]]))
    x[[4]] <- sum(GParamFlags[c("readable", "writable", "construct")])
  
  if (type == "Boolean")
    x[[5]] <- ifelse(is.null(x[[5]]), F, as.logical(x[[5]]))
  else if (type == "String")
    x[[5]] <- ifelse(is.null(x[[5]]), "", as.character(x[[5]]))
  else if (type == "Unichar")
    x[[5]] <- ifelse(is.null(x[[5]]), as.integer(0), as.integer(x[[5]]))
  else if (type %in% c("Flags", "Enum", "Param", "Boxed", "Object", "GType"))
    x[[5]] <- as.GType(x[[5]])
  else if (type %in% c("Char", "UChar")) {
    x[[5]] <- ifelse(is.null(x[[5]]), 0, as.raw(x[[5]]))
    x[[6]] <- as.raw(x[[6]])
    x[[7]] <- as.raw(x[[7]])
  } else if (type == "Int") {
    x[[5]] <- ifelse(is.null(x[[5]]), 0, as.integer(x[[5]]))
    x[[6]] <- as.integer(x[[6]])
    x[[7]] <- as.integer(x[[7]])
  } else if (type %in% c("UInt", "ULong", "Long", "UInt64", "Int64", "Float", "Double")) {
    x[[5]] <- ifelse(is.null(x[[5]]), 0, as.numeric(x[[5]]))
    x[[6]] <- as.numeric(x[[6]])
    x[[7]] <- as.numeric(x[[7]])
  } else if (type == "ValueArray")
    x[[5]] <- as.GParamSpec(x[[5]])
  else if (type == "Sexp") {
    # if there's no type, try to get it from the default value
    if (is.null(x[[5]]) && !is.null(x[[6]]))
      x[[5]] <- typeof(x[[6]])
    else if (is.null(x[[5]])) # otherwise, fallback to ANY
      x[[5]] <- "any"
    # if there's no default value, create one given the type (ANY->NULL)
    anysxp <- .RGtkCall("getNumericType", "any")
    if (is.null(x[[6]]) && x[[5]] != "any" && x[[5]] != anysxp)
      x[[6]] <- new(x[[5]])
    # if type is numeric, assume it's a type code, otherwise assume it's a type 
    # name and ask the C side to query the for the code
    if (!is.numeric(x[[5]]))
      x[[5]] <- .RGtkCall("getNumericType", x[[5]])
  }
  
  return(x)
}
  

gObjectSetData <-
function(obj, key, data = NULL)
{
        checkPtrType(obj, "GObject")
        key <- as.character(key)

        w <- .RGtkCall("S_g_object_set_data", obj, key, data, PACKAGE = "RGtk2")

        return(invisible(w))
}
gObjectGetData <-
function(obj, key)
{
        checkPtrType(obj, "GObject")
        key <- as.character(key)

        w <- .RGtkCall("S_g_object_get_data", obj, key, PACKAGE = "RGtk2")

        return(w)
}

# Methods

parentHandler <-
function(method, obj = NULL, ...)
{
  # untested stuff
  # chaining up is only allowed/makes sense from inside a GObject implementation
  stopifnot(implements(obj, "SGObject"))
  if (is.null(attr(obj, ".private")))
    stop("Parent methods should only be invoked within the instance")
  if (FALSE) { # stuff that might work some day
  parent_call <- sys.call(sys.parent(1))
  parent_frame <- parent.frame()
  formal_args <- formals(parent_call[[1]])
  if (missing(obj))
    obj <- get(names(formal_args)[1], parent_frame)
  args <- list(...)
  formal_names <- names(formal_args)[-1]
  missing_names <- formal_names[!(formal_names %in% names(args))]
  unnamed <- sapply(names(args), nchar) == 0
  names(args)[unnamed] <- missing_names[seq(along=unnamed)]
  missing_names <- missing_names[!(missing_names %in% names(args))]
  args[missing_names] <- lapply(missing_names, get, parent_frame)
  parent <- .Call("S_g_object_parent", obj, PACKAGE = "RGtk2")
  if (!is.null(parent) && is.function(try(parent[[method]], T)))
    fun <- parent[[method]]
  else { # fallback to calling a wrapper of the C virtual
    fun <- eval(substitute(gTypeGetClass(class(obj)[2])$sym, list(sym=method)))
    args <- c(obj, args)
  }      
  do.call(fun, args)
  }
  # assume looking for a function, does not make sense for fields
  #function(...) {
    # is this a function defined by a parent R class?
    parent <- .Call("S_g_object_parent", obj, PACKAGE = "RGtk2")
    if (!is.null(parent) && is.function(try(parent[[method]], T))) {
      parent[[method]](...)
    } else # fallback to calling a wrapper of the C virtual
      eval(substitute(gTypeGetClass(class(obj)[2])$sym(obj, ...), 
        list(obj=obj,sym=method)))
  #}
}

"$.<invalid>" <-
function(obj, name)
{
  stop("attempt to call '", name, "' on invalid reference '", deparse(substitute(obj)), "'", call.=FALSE)
}

"$.GObject" <- "$.RGtkObject" <-
function(x, member)
{ # try for a declared method first, else fall back to member
 result <- try(.getAutoMethodByName(x, member, parent.frame()), T)
 if (inherits(result, "try-error"))
   result <- x[[member]]
 result
}

.getAutoMemberByName <- 
function(obj, name)
{
  # if we have an SGObject, try private env (includes protected) then public
  stopifnot(implements(obj, "SGObject"))
  attrs <- attributes(obj)
  has_private <- ".private" %in% names(attrs)
  if (has_private)
    member <- try(get(name, attrs$.private), T)
  if (!has_private || inherits(member, "try-error"))
    member <- try(get(name, attrs$.public), T)
  if (is.function(member)) {
    # we need to add private env if it's not there
    if (!has_private)
      obj <- .Call("S_g_object_private", obj, PACKAGE="RGtk2")
    function(...) member(obj, ...)
  }
  else member
}

.getAutoMethodByName <-
function(obj, name, where)
{
    classes <- c(attr(obj, "interfaces"), class(obj))
    sym <- paste(tolower(substring(classes, 1, 1)), substring(classes, 2), toupper(substring(name, 1, 1)),
        substring(name,2), sep="")
    which <- sapply(sym, exists, where)

 if(!any(which))
   stop(paste("No such method", name, "for classes", paste(class(obj), collapse=", ")))

 method <- get(sym[which][1], where)
 
 function(...) method(obj, ...)
 
 #sym <- as.name(sym[which][1])
 
  # evaluate it to turn it into a function
  # and also get the correct environment
 #eval(substitute( function(...) {
 #                    sym(obj, ...)
 #                }, list(obj=obj,sym=sym)))
}

# Comparing pointers

"==.RGtkObject" <-
function(x, y) {
  identical(x, y)
}

# Fields

"$<-.GObject" <- "[[<-.GObject" <-
function(obj, member, value)
{ # first try for prop, then fall back to private env
  # this encourages the setting of properties, rather than using the back door
  result <- try(obj[member] <- value, T)
  if (inherits(result, "try-error")) {
    env <- attr(obj, ".private")
    if (is.null(env))
      stop("Cannot find '", member, "' to set in ", paste(class(obj),collapse=", "))
    protected_env <- parent.env(env)
    if (exists(member, protected_env))
      env <- protected_env
    assign(member, value, env)
  }
  obj
}

"[[.GObject" <-
function(obj, member, where = parent.frame())
{
  # check SGObject environments first, then fall back to field/property
  val <- try(.getAutoMemberByName(obj, member), TRUE)
  # check for C field (fast), then GObject prop
  if (inherits(val, "try-error"))
    val <- try(NextMethod("[[", where = where), TRUE)
  if (inherits(val, "try-error"))
    val <- try(obj$get(member), TRUE)
  if (inherits(val, "try-error"))
    stop("Cannot find '", member, "' for classes ", paste(class(obj), collapse=", "))
  val
}

"[[.RGtkObject" <-
  #
  #
  #
function(x, field, where = parent.frame())
{
  fun <- try(.getAutoElementByName(x, field, error = FALSE, where = where),
             TRUE)
  if (!inherits(fun, "try-error"))
    val <- fun(x)
  else stop("Cannot find '", field, "' for classes ",
            paste(class(x), collapse=", "))
  return(val)
}
# C field setting is not allowed
if (FALSE) {
"[[<-.RGtkObject" <-
  #
  #
  #
function(x, name, value)
{
  sym <- try(.getAutoElementByName(x, name, op = "Set", error = FALSE))
  if (!inherits(sym, "try-error"))
    val <- eval(substitute(sym(x, value), list(sym=sym)))
  else if(inherits(x, "GObject"))
   val <- x$set(name)
  else val <- sym
  return(val)
}
}
.getAutoElementByName <-
function(obj, name, op = "Get", error = TRUE, where = parent.frame())
{
  sym <- paste("S_", class(obj), op,
               toupper(substring(name, 1, 1)), substring(name, 2), sep = "")
  sym <- Find(function(x) is.loaded(x, PACKAGE = "RGtk2"), sym)
  
  if(is.null(sym)) {
    message <- paste("Could not", op, "field", name, "for classes",
                     paste(class(obj), collapse=", "))
    if(error)
      stop(message)
    else {
      v <- paste(message)
      class(v) <- "try-error"
      return(v)
    }
  }

  function(obj) .Call(sym, obj, PACKAGE = "RGtk2")
}

# This attempts to coerce an R object to an RGClosure that is understood on the C side
as.GClosure <- 
function(x)
{
  if (inherits(x, "GClosure"))
    x <- toRGClosure(x)
  else x <- as.function(x)
  class(x) <- "RGClosure"
  x
}

# This attempts to convert a C GClosure to an R closure 
# (with extra ref attribute that prevents recursion on C side)
toRGClosure <-
function(c_closure)
{
  checkPtrType(c_closure, "GClosure")
  closure <- function(...) {
    .RGtkCall("R_g_closure_invoke", c_closure, c(...), PACKAGE = "RGtk2")
  }
  attr(closure, "ref") <- c_closure
  closure
}

# virtuals for GObject
assign("GObject", c("set_property", "get_property"), .virtuals)

