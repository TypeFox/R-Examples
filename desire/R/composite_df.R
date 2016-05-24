##
## compose_df.R - compose a function and a desirability function
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

##
## Composite desirability function objects:
##
## All composite desirabilities are native R functions. They may have
## several attributes. These currently include:
##
## * composite.desc - short description of inner function
## * desire.function - original desirability
##

print.composite.desire.function <- function(x, ...) {
  message("Composite desirability: ")
  message("Inner function:")
  message("   ", attr(x, "composite.desc"))
  message("Desirability:")
  print(attr(x, "desire.function"), ...)
}

compositeDF <- function(expr, d, ...) {
  if ("composite.desire.function" %in% class(d))
    stop("Cannot recursivly composition desirabilty function.")
  sexpr <- substitute(expr)
  ## All this because we cannot 'match' the class of an expression...
  if (is.call(sexpr)) { ## Catch expressions:
    UseMethod("compositeDF", sexpr)
  } else {
    UseMethod("compositeDF", expr)
  }
}

compositeDF.call <- function(expr, d, ...) {
  expr <- substitute(expr)
  ## FIXME: If we save the parent frame, eval() cannot resolve x and
  ## we do not want to stick x into the parents env.
  ## pf <- parent.frame() # Save parent frame for evaluation of expr
  ev <- function(x, ...) {
    y <- eval(expr, envir=list(x=x))
    d(y, ...)
  }
  class(ev) <- "composite.desire.function"
  attr(ev, "composite.desc") <- paste("Expression: ", deparse(expr))
  attr(ev, "desire.function") <- d
  return(ev)
}

compositeDF.function <- function(expr, d, ...) {
  ## FIXME: merge ... of ev and ... of cdf.f:
  ev <- function(x, ...)
      d(expr(x), ...)
  class(ev) <- "composite.desire.function"
  attr(ev, "composite.desc") <- paste("Function: ", deparse(substitute(expr)), "(x)", sep="")
  attr(ev, "desire.function") <- d
  return(ev)
}

compositeDF.lm <- function(expr, d, ...) {
  ## Calculate sigma
  sigma <- summary(expr)$sigma
  ev <- function(x, ...) {
    ## Convert non data frame x arguments
    if (!is.data.frame(x)) {
      if (is.vector(x)) {
        names(x) <- pnames
        x <- as.data.frame(as.list(x))
      } else if (is.matrix(x)) {
        colnames(x) <- pnames
        x <- as.data.frame(x)
      } else {
        stop("Cannot convert argument 'x' into a data.frame object.")
      }      
    }
    y <- predict(expr, newdata=x)
    ## If this is a realistic DF, pass sd on.
    d(y, sd=sigma, ...)
  }
  ## Extract vector of names of preditor variables:
  pnames <- all.vars(formula(expr)[[3]])
  attr(ev, "composite.desc") <- paste("Linear Model: ", deparse(expr$call))
  class(ev) <- "composite.desire.function"
  attr(ev, "desire.function") <- d
  return(ev)
}
