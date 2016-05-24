### RQtLibrary objects refer to a Smoke module describing
### a bound library (e.g. Qt).

### FIXME: Think about using some S4 here

## Get available Smoke modules
qsmokes <- function() {
  .Call("qt_qsmokes", PACKAGE="qtbase")
}

qsmoke <- function(x) {
  nm <- attr(x, "name")
  smoke <- qsmokes()[[nm]]
  if (is.null(smoke))
    stop("Smoke module for library '", nm, "' not found")
  smoke
}
  
## Get classes in library
qclasses <- function(x) {
  .Call("qt_qclasses", qsmoke(x), PACKAGE="qtbase")
}

.qlibraries <- new.env(parent = emptyenv())

qlibraries <- function() as.list(.qlibraries)

## Many libraries define all of their classes within a namespace of
## the same name. We want to avoid syntax like Qanviz$Qanviz$Layer, so
## the top namespace is implied. Qt itself is of course an exception.
qlibrary <- function(lib, namespace = deparse(substitute(lib)), register = TRUE)
{
  force(namespace)
  name <- tolower(deparse(substitute(lib)))
  if (is.null(attr(lib, "name")))
    attr(lib, "name") <- name
  if (register)
    assign(attr(lib, "name"), lib, .qlibraries)
  attr(lib, "ns") <- namespace
  class(lib) <- c("RQtLibrary", "environment")
  classes <- qclasses(lib)
  if (!is.null(namespace)) { # remove the implied namespace
    prefix <- paste("^", namespace, "::", sep = "")
    classes <- grep(prefix, classes, value = TRUE)
    names(classes) <- sub(prefix, "", classes)
  } else names(classes) <- classes
  hasPrefix <- grepl("::", names(classes))
  ## take care of non-namespaced/non-internal classes first
  lapply(names(classes[!hasPrefix]), function(classAlias) {
    getClass <- function() {
      class <- qsmokeClass(lib, classes[[classAlias]])
      rm(list = classAlias, envir = lib)
      assign(classAlias, class, lib) ## cache for further use
      lockBinding(classAlias, lib)
      class
    }
    makeActiveBinding(classAlias, getClass, lib)
  })
  ## now we need a separate environment for each namespace
  ns <- sub("(.*?)::.*", "\\1", names(classes)[hasPrefix])
  if (length(ns) && !is.null(namespace))
    ns <- paste(namespace, ns, sep = "::")
  ns <- setdiff(ns, classes)
  for (nsi in ns) {
    env <- new.env()
    attributes(env) <- attributes(lib)
    assign(sub(paste("^", namespace, "::", sep = ""), "", nsi),
           qlibrary(env, nsi, FALSE), lib)
  }
  lib
}

## Usually, one just uses '$' to lookup classes; this will map a full
## class name through nested namespaces and internal classes.
qclassForName <- function(name, lib) {
  classes <- vector("list", length(name))
  if (missing(lib)) {
    libs <- qlibraries()
    notfound <- rep(TRUE, length(classes))
    for(lib in libs) { # simple linear search probably OK for now
      classes[notfound] <- qclassForName(name[notfound], lib)
      notfound[notfound] <- unlist(lapply(classes[notfound], is.null))
    }
    return(classes)
  }
  if (is.null(lib))
    return(classes)
  .mget <- function(x) mget(x, lib, ifnotfound = list(NULL))
  if (is(lib, "RQtClass")) {
    namespace <- attr(lib, "name")
    lib <- attr(lib, "env")
  }
  else namespace <- attr(lib, "ns")
  name <- sub(paste(namespace, "::", sep = ""), "", name)
  hasPrefix <- grepl("::", name)
  classes[!hasPrefix] <- .mget(name[!hasPrefix])
  ns <- factor(sub("^([^:])*::.*", "\\1", name[hasPrefix]))
  if (length(ns)) # unsplit() broken for zero-length factors
    classes[hasPrefix] <- unsplit(mapply(qclassForName,
                                         split(name[hasPrefix], ns),
                                         .mget(levels(ns))),
                                  ns)
  classes
}

print.RQtLibrary <- function(x, ...) {
  cat("Module '", attr(x, "name"), "' with ", length(ls(x)),
      " top-level classes\n", sep = "")
}

### Module object for Qt library.
Qt <- new.env()

