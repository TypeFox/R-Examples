#' Define a new class
#' 
#' This is an experimental implementation of reference classes. Use \code{\link{defineRefClass}} or \code{\link{retList}} instead. \code{defineClass} has side effects. The constructor is the return value of \code{defineClass}.
#'  
#' @param name character name of the class
#' @param expr expression
#' @param contains character name of class from which to inherit
#' 
#' @details 
#' \code{defineClass} creates a S4-Class which can be used for standard S4 method dispatch. It will also set the method 'initialize' which need not to be changed. If you want to have some operations carried out on initialization use a function definition named \code{init} as part of \code{expr}. The return value from \code{defineClass} is the constructor function. It has the argument \code{...} which will be passed to \code{init}.
#' 
#' All classes defined with \code{defineClass} inherit from class "aoos" which is a S4-class containing an environment. In that environment \code{expr} is evaluated; for inheritance, all \code{expr} from all parents will be evaluated first.
#' 
#' Everything in \code{expr} will be part of the new class definition. A leading dot in a name will be interpreted as private. You can use \code{public} and \code{private} to declare private and public members explicitly. If \code{x} in a call to \code{public} is a function it will be a public member function (method). For any other class the return value of \code{public} is a get and set method. If called without argument it will get the value, if called with argument it will set the value. You can define a validity function which will be called whenever the set method is called. Objects which inherit from class \code{environment} can be accessed directly, i.e. not via get/set methods. If you want to access fields without get/set methods, you can use the class \code{\link{Accessor-class}}.
#' 
#' @seealso \code{\link{Accessor-class}}, \code{\link{Binary-class}}, \code{\link{Show-class}}
#' 
#' @rdname defineClass
#' @export
#' @examples
#' test <- defineClass("test", {
#'   x <- "Working ..."
#'   .y <- 0
#'   doSomething <- public(function() {
#'     self$.y <- .y + 1
#'     cat(x(), "\n")
#'     invisible(self)
#'   })
#' })
#' instance <- test()
#' \dontrun{
#' instance$.y # error
#' }
#' instance$doSomething()$doSomething()
#' instance$x()
#' instance$x(2)
#' instance$x()
#'
#' # Example for reference classes as field
#' MoreTesting <- defineClass("MoreTesting", {
#'   refObj <- test()
#' })
#' instance <- MoreTesting()
#' instance$refObj$x()
defineClass <- function(name, expr, contains = NULL) {
  
  expr <- substitute(expr)
  parentEnv <- parent.frame()
  aoosClass <- findAoosClasses(contains)
  
  getMember <- function() {
    e <- setEnvironment(aoosClass, parentEnv)
    eval(expr, envir = e)
    arrangeEnvironment(e)
  }
  
  const <- function(...) {
    do.call("new", c(list(Class = name), ...))
  }
  
  setClass(name, where = parentEnv, 
           contains = c(contains, if(!length(aoosClass)) "aoos"))
  
  setMethod("initialize", name,
            function(.Object, ...) {
              .Object@.xData <- getMember()
              parent.env(.Object)$self <- .Object
              parent.env(.Object)$.self <- .Object
              init(.Object, ...)
            }, where = parentEnv)
  
  invisible(const)
}

findAoosClasses <- function(contains) {
  ind <- sapply(contains, function(class) {
    any(names(getClassDef(class)@contains) == "aoos")
  })
  contains[ind]
}

setEnvironment <- function(aoosClass, parentEnv) {
  if(!length(aoosClass)) {
    new.env(parent = parentEnv)
  } else {
    object <- new(aoosClass)
    parent.env(object)
  }
}

arrangeEnvironment <- function(e) {
  handleSpecialNames(e)
  handlePublicNames(e)
  f <- makePublicRepresentation(e)
  parent.env(f) <- e
  f
}

handleSpecialNames <- function(e) {
  if(exists("init", envir = e, inherits = FALSE)) 
    assign("init", private(get("init", envir = e)), envir = e)
}

handlePublicNames <- function(e) {
  publicNames <- ls(envir = e)
  lapply(publicNames, function(n) assign(n, public(get(n, envir = e)), envir = e))
}

makePublicRepresentation <- function(e) {
  allMember <- as.list(e, all.names = TRUE)
  publicMemberInd <- sapply(allMember, function(obj) inherits(obj, "public"))
  publicMember <- lapply(allMember[publicMemberInd], getPublicRepresentation)
  list2env(publicMember)
}

init <- function(object, ...) {
  
  shouldInitBeCalled <- function() {
    initMethodExists <- exists("init", envir = parent.env(object), inherits = FALSE)
    #hasArguments <- as.logical(length(list(...)))
    initMethodExists
  }
  
  if(shouldInitBeCalled()) {
    parent.env(object)$init(...)
  }
  
  object
}
