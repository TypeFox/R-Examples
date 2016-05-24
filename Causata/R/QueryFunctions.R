# functions that operate on Query and FocalPointQuery
#
# author: davidb <support@causata.com>
library(stringr)

#
# WithVariables
#
WithVariables <- function(...) {
  this <- list(variables=as.character(c(list(...), recursive=TRUE)))
  class(this) <- "WithVariables"
  this
}
is.WithVariables <- function(this) inherits(this, "WithVariables")
Variables.WithVariables <- function(this) {
  this$variables
}


#
# Where
#
# Creates a clause than can be added to a Query or FocalPointQuery
#
# Example uses:
# Where("`variable-name-no-escaping` > 123")
#   Creates a where clause: WHERE `variable-name-no-escaping` > 123
# Where("name.will.be.escaped", "> 123")
#   Creates a where clause: WHERE variable.`name.will.be.escaped`>123
# Where("name.will.be.escaped", GreaterThan(123)
#   Creates a where clause: WHERE `name.will.be.escaped`>123
Where <- function(variable, operation=NULL) {
  if (is.null(operation)) {
    # single argument, assume it's a clause
    this <- list(clause=variable)
  } else {
    # two arguments
    filter <- operation
    if (is.RawOperator(filter)) {
      this <- list(clause=paste(c("variable.", Backtick(variable), " ", filter$name, " ", filter$operand), sep="", collapse=""))
    } else {
      this <- list(clause=paste(c("variable.", Backtick(variable), " ", filter), collapse=""))
    }
  }
  
  class(this) <- "Where"
  return(this)
}
is.Where <- function(this) inherits(this, "Where")

#
# Limit
#
Limit.numeric <- function(this, ...) {
  obj <- list(limit=this)
  class(obj) <- "Limit"
  return(obj)
}
is.Limit <- function(this) inherits(this, "Limit")


#
# Operators
#
EscapeQuotes <- function(value) str_replace_all(value, "'", "\\\\'")
QuoteString <- function(value) paste("'", EscapeQuotes(value), "'", sep="")
QuotedIfCharacter <- function(value) {
  if (class(value) == "character") {
    QuoteString(value)
  } else {
    value
  }
}

RawOperator <- function(name, operand) {
  this <- list(name=name, operand=operand)
  class(this) <- "RawOperator"
  this
}
is.RawOperator <- function(this) inherits(this, "RawOperator")

# TODO: handle dates and convert them to ISO format (which ISO format, I don't know)
Operator <- function(symbol, rhs) {
  RawOperator(symbol, as.character(QuotedIfCharacter(rhs)))
}

NumericOperator <- function(symbol, rhs) {
  stopifnot(is.numeric(rhs))
  Operator(symbol, rhs)
}

TextOperator <- function(symbol, rhs) {
  stopifnot(is.character(rhs))
  Operator(symbol, rhs)
}

EqualTo                <- function(value) Operator("=", value)
NotEqualTo             <- function(value) Operator("<>", value)
GreaterThan            <- function(value) NumericOperator(">", value)
GreaterThanOrEqualTo   <- function(value) NumericOperator(">=", value)
LessThan               <- function(value) NumericOperator("<", value)
LessThanOrEqualTo      <- function(value) NumericOperator("<=", value)
Like                   <- function(value) TextOperator("LIKE", value)
BeginningWith          <- function(value) TextOperator("BEGINNING WITH", value)
Between                <- function(value) {
  stop("BETWEEN not yet supported")
}
In <- function(...) {
  values <- paste(sapply(as.character(c(list(...), recursive=TRUE)), QuoteString), collapse=",")
  operator.rhs <- paste("(", values, ")", sep="")
  RawOperator("IN", operator.rhs)
}


BacktickCollapse <- function(variable.names){
  # converts a list of variable names to SQL format with backticks, returns a single string
  return(paste(Backtick(variable.names), collapse=","))
}


# private
# adds back ticks to the items passed in.
# whatever is passed is converted to a vector and returned as a vector with backticked items
# > Backtick(c("a", "b"))
# [1] "`a`" "`b`"
# > Backtick(list("a", "b"))
# [1] "`a`" "`b`"
# > Backtick("a", "b")
# [1] "`a`" "`b`"
#
Backtick <- function(...) {
  backticker <- function(value) paste("`", value, "`", sep="")
  as.vector(sapply(as.character(c(list(...), recursive=TRUE)), backticker))
}
