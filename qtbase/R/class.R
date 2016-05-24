### Classes are functions (constructors) with static methods in an environment

## invoke a static method
"$.RQtClass" <- function(x, name) {
  attr(x, "env")[[name]]
}
"[[.RQtClass" <- function(x, name) {
  attr(x, "env")[[name]]
}


names.RQtClass <- function(x) {
  ls(attr(x, "env"))
}

qmethods <- function(x) {
  stopifnot(is(x, "RQtClass"))
  methods <- .Call("qt_qmethods", x, PACKAGE="qtbase")
  names(methods) <- c("name", "return", "signature", "protected", "static",
                      "constructor")
  df <- as.data.frame(methods, stringsAsFactors=FALSE)
  df[!duplicated(df$signature),]
}

qenums <- function(x) {
  stopifnot(is(x, "RQtClass"))
  .Call("qt_qenums", x, PACKAGE="qtbase")
}

qparentClasses <- function(x) {
  stopifnot(is(x, "RQtClass"))
  .Call("qt_qparentClasses", x, PACKAGE="qtbase")
}

## We enforce single inheritance in the implementation of
## RQtUserClass, but we abstract that here.
qparents <- function(x, ...) UseMethod("qparents")
qparents.RQtSmokeClass <- function(x) attr(x, "parents")
qparents.RQtUserClass <- function(x) {
  parent <- attr(x, "parent")
  structure(list(parent), names = attr(parent, "name"))
}

print.RQtClass <- function(x, ...) {
  methods <- qmethods(x)
  public <- methods[!methods$protected,]
  cat("Class '", attr(x, "name"), "' with ", nrow(public), " public methods\n",
      sep = "")
}

## Smoke classes are populated with their entire hierarchy, as
## environment inheritance is single, while C++ inheritance is
## multiple.

## obtain a class object from a smoke module and a name
qsmokeClass <- function(x, name, internals = character()) {
  env <- new.env(parent = emptyenv())
  ## we have to be extra careful about the enclosing environments for
  ## these functions; otherwise, we could end up e.g. serializing the
  ## big library 'x' -- could do this more idiomatically by defining
  ## the functions outside the scope of this one.. but lets just be explicit.
  constructor <- function(...) qinvokeStatic(cl, basename, ...)
  constructorEnv <- new.env(parent = getNamespace("qtbase"))
  environment(constructor) <- constructorEnv
  constructorEnv$basename <- gsub(".*::", "", name)
  cl <- structure(constructor, name = name,
                  env = env, module = attr(x, "name"),
                  class = c("RQtSmokeClass", "RQtClass", "function"))
  constructorEnv$cl <- cl
  attr(cl, "parents") <- qclassForName(qparentClasses(cl))
  methods <- qmethods(cl)
  methods <- methods[!duplicated(methods$name) & methods$static &
                         !methods$protected,]
  lapply(methods$name, function(name) {
    fun <- structure(function(...) qinvokeStatic(cl, name, ...), static = TRUE)
    environment(fun) <- list2env(list(cl = cl, name = name),
                                 parent = getNamespace("qtbase"))
    assign(name, fun, env)
  })
  enums <- qenums(cl)
  for (enum in names(enums))
    assign(enum, structure(enums[enum], class = "QtEnum"), env)
  internals <- grep(paste("^", name, "::", sep = ""), qclasses(x), value = TRUE)
  for (internal in internals)
    assign(gsub(".*::", "", internal), qsmokeClass(x, internal), env)
  lockEnvironment(env, TRUE)
  cl
}

## Rewrites the constructor to handle base initialization
normConstructor <- function(x, parent) {
  if (!is.function(x))
    stop("constructor must be a function")

  ## Find base constructor call (the first one)
  b <- body(x)
  first <- which(sapply(b, identical, as.name("{"))) + 1
  expr <- NULL
  if (!length(first)) # no braces
    expr <- b
  else if (first <= length(b))
    expr <- b[[first]]
  if (is.call(expr) && identical(expr[[1]], as.name("super"))) {
    baseCall <- expr
    baseCall[[1]] <- parent
    if (length(first))
      body(x)[[first]] <- NULL
  } else baseCall <- as.call(list(parent)) # make default
  
  ## Stick the call into its own function, so we can enclose it
  baseConstructor <- as.function(c(formals(x), baseCall))
  environment(baseConstructor) <- environment(x)
  
  ## Now generate a new function that:
  ## - executes the base constructor
  ## - casts the base instance down to this class
  ## - encloses the constructor in the instance env (enclosed by original env)
  ## - invokes the constructor
  argNames <- lapply(names(formals(x)), as.name)
  wrapperBody <- substitute({
    base <- callBaseConstructor
    this <- qcast(base, sys.function())
    x_enclosed <- qenclose(this, x)
    callConstructor
    this
  }, list(callBaseConstructor = as.call(c(baseConstructor, argNames)),
          callConstructor = as.call(c(quote(x_enclosed), argNames)), x = x))
  if (!length(first)) # no need to call constructor
    wrapperBody[[5]] <- NULL
  fun <- as.function(c(formals(x), wrapperBody))
  environment(fun) <- getNamespace("qtbase")
  fun
}

qsetClass <- function(name, parent, constructor = function(...) parent(...),
                      where = topenv(parent.frame()))
{
  ## mangle the class name to prevent conflicts
  module <- getPackageName(where)
  prefixedName <- paste("R", module, name, sep = "::")
  if (exists(name, where))
    warning("Symbol '", name, "' already exists in '", module, "'")
  ## get our real constructor
  constructor <- normConstructor(constructor, parent)
  ### FIXME: May want to support reregistration of classes. This requires:
  ### 1) chaining up at the C++ Class level, rather than at instanceEnv
  ### 2) reducing the 'parent' attribute to a light-weight reference
  ### 3) add ability to unregister classes from cache
  ### NOTE: Not sure if this is a good idea, since it breaks instances
  parentEnv <- attr(parent, "instanceEnv")
  if (is.null(parentEnv))
    parentEnv <- emptyenv() # a smoke class, no instance symbols
  instanceEnv <- new.env(parent = parentEnv)
  env <- attr(parent, "env") # do not support user static methods yet
  metadata <- new.env(parent = emptyenv())
  metadata$properties <- new.env(parent = emptyenv())
  cl <- structure(constructor, module = module, name = prefixedName,
                  parent = parent, env = env, instanceEnv = instanceEnv,
                  metadata = metadata,
                  class = c("RQtUserClass", "RQtClass", "function"))
  qinitClass(cl)
  assign(name, cl, where)
  cl
}

## ensures that the class is (re)initialized
qinitClass <- function(x) {
  .Call("qt_qinitClass", x, PACKAGE="qtbase")
}

qcast <- function(x, class) {
  .Call("qt_qcast", x, class, PACKAGE="qtbase")
}

qenclose <- function(x, fun) {
  .Call("qt_qenclose", x, fun, PACKAGE="qtbase")
}

qsetMethod <- function(name, class, FUN,
                       access = c("public", "protected", "private"))
{
  attr(FUN, "access") <- match.arg(access)
  assign(name, FUN, attr(class, "instanceEnv"))
  name
}

qhasMethod <- function(name, class) {
  exists(name, attr(class, "instanceEnv"))
}

## Integration with the Qt Meta Object Compiler (MOC)

## The basic idea: define methods in R that are described by
## QMetaObject. This allows R to define signals and slots (and
## properties, enums, etc). The utility of signals is obvious. Slots
## could be exposed as e.g. DBus services. The main downside is that
## providing an external interface requires us to specify the types
## using C++ nomenclature.

## The methods will belong to the R class, as usual. We will compile a
## QMetaObject and provide it via the QObject::metaObject() virtual
## method. Then we will catch invocations via the QObject::qt_metacall
## virtual. All methods could be forwarded to R, but we might
## short-circuit signal emissions (and call QMetaObject::activate).

qsetSlot <- function(signature, class, FUN,
                     access = c("public", "protected", "private"))
{
  access <- match.arg(access)
  method <- qmetaMethod(signature, access, names(formals(FUN)))
  qsetMethod(method$name, class, FUN, access)
  qmetadata(class)$slots[[signature]] <- method
  method$name
}

## Signals are essentially implemented by QMetaObject::activate(). We
## could have the signal method call this directly, but for
## convenience we instead call the corresponding QMetaMethod. This is
## caught by the qt_metacall override which then calls
## QMetaObject::activate().

qsetSignal <- function(signature, class,
                       access = c("public", "protected", "private"))
{
  access <- match.arg(access)
  method <- qmetaMethod(signature, access)
  meta <- qmetaObject(class)
  index <- meta$methodCount()
  qsetMethod(method$name, class,
             function(...) .Call("qt_qmetaInvoke", this, index, list(...),
                                 PACKAGE="qtbase"),
             access)
  ## set this last so that the compiled metadata does not end up in a
  ## lazy-loaded package
  qmetadata(class)$signals[[signature]] <- method
  method$name
}

qsetProperty <- function(name, class, type = NULL,
                         read = function() this[[.name]],
                         write = function(val) this[[.name]] <- val,
                         ##reset = NULL,
                         notify = NULL, 
                         constant = FALSE, final = FALSE,
                         ##designable = TRUE, scriptable = TRUE,
                         stored = TRUE, user = FALSE)
{
  ## FIXME: obviously have to do better job of checking arguments here
  ## FIXME: do we want to map R types to C++ types for 'type'?
  ##        - if so, we also need to map method signatures
  if (missing(name) || !is.character(name))
    stop("'name' is required, as character vector")
  if (missing(class))
    stop("'class' is required")
  if (is.null(type)) {
    if (!missing(constant) || !missing(final) || !missing(stored) ||
        !missing(user))
      stop("Arguments 'constant', 'final', 'stored', and 'user' are",
           " ignored if 'type' is NULL")
  } else if (!is.character(type))
    stop("'type' should be NULL or a character vector")
  if (!is.null(notify)) {
    notify <- qresolveSignature(class, notify, "signal")
    writeArg <- formals(write)
    writeExpr <- call("function", writeArg, body(write))
    writeBody <- call("{",
                      as.call(list(writeExpr, as.name(names(writeArg)))),
                      call(sub("\\(.*", "", notify)))
    write <- as.function(c(writeArg, writeBody), environment(write))
  }
  .name <- paste(".", name, sep = "")
  prop <- list(name = name, type = type, read = read,
               write = write, notify = notify, constant = constant,
               final = final, stored = stored, user = user)
  qmetadata(class)$properties[[name]] <- prop
  name
}

qsetRefClass <- function(Class, where = topenv(parent.frame()), ...) {
  if (!is(Class, "RQtClass"))
    stop("'Class' must be an RQtClass, e.g., Qt$QWidget")
  parents <- qparents(Class)
  for (parent in parents)
    if (!isClass(attr(parent, "name")))
      qsetRefClass(parent, where = where)
  getPropertyNames <- function(x) rownames(qproperties(x))
  propertyNames <- getPropertyNames(Class)
  propertyNames <- setdiff(propertyNames,
                           unlist(lapply(parents, getPropertyNames)))
  fields <- sapply(propertyNames, function(propertyName) {
    eval(substitute(function(value) {
      if (missing(value))
        .ref$propertyName
      else .ref$propertyName <- value
    }, list(propertyName = as.name(propertyName))))
  })
  className <- attr(Class, "name")
  getMethodNames <- function(x) {
    methodInfo <- qmethods(x)
    methodInfo$name[!methodInfo$static & !methodInfo$constructor]
  }
  methodNames <- getMethodNames(Class)
  methodNames <- setdiff(methodNames, unlist(lapply(parents, getMethodNames)))
  methods <- sapply(methodNames, function(methodName) {
    eval(substitute(function(...) {
      qinvoke(.ref, methodName, ...)
    }, list(methodName = methodName)))
  })
  methods <- c(methods,
               initialize = eval(substitute(function(...) {
                 Class(...)
               }, list(Class = Class))))
  
  setRefClass(className, fields = c(fields, .ref = "RQtObject"),
              methods = methods, contains = sapply(parents, attr, "name"),
              where = where, ...)
}
