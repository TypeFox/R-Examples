
############## General ######################
#' Prints information of a \link{ml.col.def-class} object.
#'
#' @param x an \link{ml.col.def-class} object
#' @export
setMethod("print", signature(x="ml.col.def"),
          function (x) {
            cat(paste("Column definition: ", x@.expr,
                      "\nMarkLogic element/property name:", x@.org_name, "\nSource document format:", x@.format,
                      "\nSource document namespace:", x@.xmlns,
                      "\nUse this to define new columns on a ml.data.frame using the $ operator. To select a subset of a ml.data.frame, use bldf[] notation. "),"\n")
          }
)
#' Prints information of a \link{ml.col.def-class}
#'
#' @param object an \link{ml.col.def-class} object
#' @export
setMethod("show", signature(object="ml.col.def"),
          function (object) {
            cat(paste("Column definition:", object@.expr,
                      " \nMarkLogic element/property name:", object@.org_name, "\nSource document format:", object@.format,
                      "\nSource document namespace:", object@.xmlns
                      ,"\nUse this to define new columns on a ml.data.frame using the $ operator. To select a subset of a ml.data.frame, use bldf[] notation. "),"\n")
          }
)
#' Check if an object is of type \link{ml.col.def-class}
#'
#' This function checks if the input is of type \link{ml.col.def-class}.
#'
#' @param x The input can be of any type.
#' @return True if it is a \link{ml.col.def-class}. False otherwise.
#' @export
is.ml.col.def <-
  function(x) {
    return(inherits(x, "ml.col.def"))
  }
# internal used function
as.ml.col.def <- function(x) {
  if(inherits(x,"ml.data.frame")) {
    return(paste('"rfmlResult[\'',x@.col.name[1],'\']"',sep=''))
  } else if(inherits(x,"ml.col.def")) {
    return(x@.expr)
  } else {
    if(is.character(x)) {
      return(paste("'",x,"'",sep=''))
    } else {
      return(as.character(x))
    }
  }
}

################ Arithmetic operators ############################
#' Arithmetic Operators
#'
#' @param e1,e2 numeric vectors or string or \link{ml.col.def-class} object.
#' @name arith
NULL

#' @rdname arith
#' @export
setMethod("Arith", signature(e1="ml.col.def",e2="ml.col.def"), function(e1, e2) {
  if(e1@.parent@.name!=e2@.parent@.name) {
    stop("Cannot combine two columns from different ml.data.frames.")
  }
  # need to check the data types of
  if (e1@.data_type == e2@.data_type) {
    dataType = e1@.data_type
  } else {
    # we only use string and number, so if they are not same we fall back to string
    dataType = "string"
  }
  return(new(Class="ml.col.def",
             .expr=paste('(',as.ml.col.def(eval(e1)),.Generic,as.ml.col.def(eval(e2)),')',sep=''),
             .data_type=dataType,
             .parent=e1@.parent,
             .type="expr",
             .aggType=aggType(e1,e2)));
})

#' @rdname arith
#' @export
setMethod("Arith", signature(e1="ml.col.def", e2="ANY"), function(e1, e2) {
  if (!(is.numeric(e2)) || e1@.data_type == "string") {
    dataType <- "string"
  } else {
    dataType <- "number"
  }
  return(new(Class="ml.col.def",.expr=paste('(', as.ml.col.def(eval(e1)),.Generic,as.ml.col.def(eval(e2)),')',sep=''),.data_type=dataType,.parent=e1@.parent,.type="expr",.aggType=aggType(e1,e2)));
})

#' @rdname arith
#' @export
setMethod("Arith", signature(e1="ANY", e2="ml.col.def"), function(e1, e2) {
  if (!(is.numeric(e1)) || e2@.data_type == "string") {
    dataType <- "string"
  } else {
    dataType <- "number"
  }
  return(new(Class="ml.col.def",.expr=paste('(', as.ml.col.def(eval(e1)),.Generic,as.ml.col.def(eval(e2)),')',sep=''),.data_type=dataType,.parent=e2@.parent,.type="expr",.aggType=aggType(e1,e2)));
})

################ Comparison operators ############################
#' Relational Operators
#'
#' Relational operators used for field level filtering of a \link{ml.data.frame} object.
#'
#' @param e1 an \link{ml.col.def-class} object.
#' @param e2 any object
#' @export
setMethod("Compare", signature(e1="ml.col.def", e2="ANY"), function(e1, e2) {
  mlDf <- e1@.parent
  i <- which(mlDf@.col.name %in% e1@.name)
  # need to handle the comparsion operator...
  fieldQuery <- paste('{"', mlDf@.col.org_name[i],
                     '":{"value":"',e2,
                     '","operator":"',.Generic,
                     '","orgPath":"',mlDf@.col.org_xpath[i],
                     '","orgFormat":"',mlDf@.col.format[i],
                     '","xmlns":"',mlDf@.col.xmlns[i], '"}}',sep='')
  return(new(Class="ml.col.def",
             .expr=fieldQuery,
             .data_type="string",
             .parent=e1@.parent,
             .type="logical",.aggType="none"));
})

################ Scalar functions ############################
# Math functions, the following is included in Math:
# abs, sign, sqrt, floor, ceiling, trunc, round, signif,
# exp, log, expm1, log1p,cos, sin, tan,cospi, sinpi, tanpi,
# acos, asin, atan, cosh, sinh, tanh, acosh, asinh, atanh,
# lgamma, gamma, digamma, trigamma, cumsum, cumprod, cummax, cummin,
#
# We only currently support part of these, se below.

#' Miscellaneous Mathematical Functions
#'
#' Mathematical functions that can be used on \link{ml.data.frame} fields. The function is applied
#' when the result is returned to the client.
#' Only abs, acos, asin, atan, ceiling, cos, cosh, exp, floor, log, log10 ,
#' tan, tanh, sqrt, sin, sinh and trunc is currently supported.
#'
#' @param x an \link{ml.col.def-class} object.
#' @export
setMethod("Math",signature(x='ml.col.def'),function (x) {

  func <- switch(.Generic,
         abs='fn.abs',
         acos='math.acos',
         asin='math.asin',
         atan='math.atan',
         ceiling='math.ceil',
         cos='math.cos',
         cosh='math.cosh',
         exp='math.exp',
         floor='math.floor',
         log='math.log',
         log10='math.log10',
         tan='math.tan',
         tanh='math.tanh',
         sqrt='math.sqrt',
         sin='math.sin',
         sinh='math.sinh',
         trunc='math.trunc',
        stop(paste(.Generic, " not supported"))
  )

  return(new(Class="ml.col.def",.expr=paste(func, '(', as.ml.col.def(eval(x)),')',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
});

#' Cotangent
#'
#' Returns the cotangent of x.
#'
#' @param x an \link{ml.col.def-class} object.
#' @return The cotangent of x.
#' @export
cot <- function (x) {

  func <- "math.cot"
  return(new(Class="ml.col.def",.expr=paste(func, '(', as.ml.col.def(eval(x)),')',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
}

#' Degrees
#'
#' Returns numeric expression converted from radians to degrees. The function is applied
#' when the result is returned to the client.
#'
#' @param x an \link{ml.col.def-class} object.
#' @return numeric expression converted from radians to degrees.
#' @export
degrees <- function (x) {

  func <- "math.degrees"
  return(new(Class="ml.col.def",.expr=paste(func, '(', as.ml.col.def(eval(x)),')',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
}

#' Radians
#'
#' Returns numeric expression converted from degrees to radians. The function is applied
#' when the result is returned to the client.
#'
#' @param x an \link{ml.col.def-class} object.
#' @return numeric expression converted from degrees to radians.
#' @export
radians <- function (x) {

  func <- "math.radians"
  return(new(Class="ml.col.def",.expr=paste(func, '(', as.ml.col.def(eval(x)),')',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
};


################ Casting operators ############################
#' Cast a \link{ml.col.def-class} expresion to numeric.
#'
#' This function will add a function to cast the expresion of the \link{ml.col.def-class} to
#' a numeric value. The cast will be done when the result is returned to the client.
#'
#' @param  x an \link{ml.col.def-class} object.
#' @export
setMethod('as.numeric',signature(x="ml.col.def"),function (x) {
  #checkLogical(F,x);
  return(new(Class="ml.col.def",.exor=paste('Number(',as.ml.col.def(eval(x)),')',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
})

#' Cast a \link{ml.col.def-class} expresion to string
#'
#' This function will add a function to cast the expresion of the \link{ml.col.def-class} to
#' a string value. The cast will be done when the result is returned to the client.
#'
#' @param  x an \link{ml.col.def-class} object
#' @export
#' @export
setMethod('as.character',signature(x="ml.col.def"),function (x) {
  #checkLogical(F,x);
  return(new(Class="ml.col.def",.expr=paste('String(',as.ml.col.def(x),')',sep=''),.parent=x@.parent,.data_type="string",.type="expr",.aggType=aggType(x)));
})

#' Cast a \link{ml.col.def-class} expresion to integer
#'
#' This function will add a function to cast the expresion of the \link{ml.col.def-class} to
#' a integer value. The cast will be done when the result is returned to the client.
#'
#' @param  x an \link{ml.col.def-class} object
#' @export
setMethod('as.integer',signature(x="ml.col.def"),function (x) {

  #checkLogical(F,x);
  return(new(Class="ml.col.def",.expr=paste('parseInt(',as.ml.col.def(eval(x)),' )',sep=''),.parent=x@.parent,.data_type="number",.type="expr",.aggType=aggType(x)));
})

################ Utilities ############################

aggType <- function(...) {

  args = list(...);

  for(i in 1:length(args)) {
    if(inherits(args[[i]],'ml.col.def')) {
      if(args[[i]]@.aggType=="simple") {
        return("simple")
      }
    }
  }
  return("none");
}
