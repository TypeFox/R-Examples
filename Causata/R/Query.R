# Author: David Barker <support@causata.com>

Query <- function() {
  this <- list()
  class(this) <- "Query"
  this$variables <- NULL
  this$filters <- vector(mode="character")
  this$limit <- NULL
  return(this)
}
is.Query <- function(this) inherits(this, "Query")

#
# define generic functions
#
#as.character <- function(this, ...) {
#  UseMethod("as.character", this)
#}
Variables <- function(this, ...) {
  UseMethod("Variables")
}
Limit <- function(this, ...) {
  UseMethod("Limit")
}
`Limit<-` <- function(this, value) {
  UseMethod("Limit<-")
}
`Variables<-` <- function(this, value) {
  UseMethod("Variables<-")
}
#Ops <- function(this, ...) {
#  UseMethod("Ops")
#}

as.character.Query <- function(x, ...) {
  column.string <- if (length(x$variables)) paste(Backtick(x$variables), collapse=",") else "*"
  sql <- paste("SELECT", column.string, "FROM Customers variable", collapse=" ")
  
  if (length(x$filters)) {
    where.conditions <- paste(x$filters, collapse=" AND ")
    sql <- paste(sql, "WHERE", where.conditions, collapse=" ")
  }
  
  if (length(x$limit)) {
    # use sprintf with limit to prevent scientific notation (10k = 1e4) which breaks SQL
    sql <- paste(sql, "LIMIT", sprintf("%i",x$limit), collapse=" ")
  }
  
  return(sql)
}

Limit.Query <- function(this, ...) {
  this$limit
}

`Limit<-.Query` <- function(this, value) {
  this$limit <- value
  this
}

Variables.Query <- function(this, ...) {
  this$variables
}

`Variables<-.Query` <- function(this, value) {
  this$variables <- value
  this
}


# private
AddToQuery <- function(query, addition) {
  AddVariables <- function() {
    query$variables <- unique(c(query$variables, Variables(addition)))
    return(query)
  }
  AddVariableFilter <- function() {
    query$filters <- c(query$filters, addition$clause)
    return(query)
  }
  SetLimit <- function() {
    query$limit <- addition$limit
    return(query)
  }
  stopifnot(is.Query(query))
  if (is.WithVariables(addition)) {
    AddVariables()
  } else if (is.Where(addition)) {
    AddVariableFilter()
  } else if (is.Limit(addition)) {
    SetLimit()
  } else {
    stop(paste("Cannot add ", class(addition), " to Query"))
  }
}

Ops.Query <- function(e1, e2) {
  lhs <- e1 # copy
  rhs <- e2
  if (nargs() == 1) stop("Unary ", .Generic, " not defined for Query")
  if (.Generic != "+") stop("Only + operator is defined for Query")
  
  if (is.Query(lhs)) {
    AddToQuery(lhs, rhs)
  } else if (is.WithVariables(rhs)) {
    AddToQuery(rhs, lhs)
  } else if (is.Where(rhs)) {
    AddToQuery(rhs, lhs)
  } else if (is.Limit(rhs)) {
    AddToQuery(rhs, lhs)
  } else {
    NextMethod(.Generic)
  }
}
