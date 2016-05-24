## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## WORKAROUND for namespace clashes between R.oo and rJava
## where some R.oo S3 methods for Object and Exception
## override the intended ones for rJava objects with
## class attributes containing these classes as well.
##
## See https://github.com/s-u/rJava/issues/60
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
.fixMethodS3 <- function(generic, class, expr=NULL, envir=parent.frame()) {
  method <- sprintf("%s.%s", generic, class)
  expr <- substitute(expr)

  f <- get(method, mode="function", envir=getNamespace("R.oo"), inherits=TRUE)
  if (is.null(expr)) {
    x <- as.symbol(names(formals(f)[1]))
    expr <- substitute(
      if(!.isRoo(x)) return(NextMethod())
    , list(x=x))
  }
  
  body(f) <- substitute({
    a
    b
  }, list(a=expr, b=body(f)))
  assign(method, f, envir=envir, inherits=TRUE)
 
  invisible(f)
} ## .fixMethodS3()

.isRoo <- function(x) is.environment(attr(x, ".env"))

.fixMethodS3("names", "Object")
.fixMethodS3("$", "Object")
.fixMethodS3("[[", "Object")
.fixMethodS3("print", "Object")
.fixMethodS3("print", "Exception")
