# Author: David Barker <support@causata.com>

FocalPointQuery <- function(
  focalpoint.event,
  cardinality=if (length(event.attribute)) { "using.all.values" } else { "using.all.events" },
  event.attribute=NULL) {
  # Legal cardinalities,
  # when just specifying an event, not an event attribute:
  #   using.all.events, using.oldest.event, using.newest.event
  # when specifying an event attribute, the focalpoint points(s) are taken from the values of that attribute.
  # Therefore, the type of the event attribute must be a DATE
  #   using.all.values, using.earliest.value, using.most.recent.value
  #
  # For events, these cardinalities map to SQL like this:
  # using.oldest.event -> AND IS_LAST(variable.timestamp)
  # using.newest.event -> AND IS_FIRST(variable.timestamp)
  #
  #
  this <- list(
    cardinality=as.character( cardinality ),
    focalpoint.event=as.character( focalpoint.event ),
    event.attribute=if (is.null( event.attribute )) NULL else as.character( event.attribute ),
    variables=NULL,
    where.clauses=vector(mode="character"),
    extra.events=vector(mode="character"),
    limit=NULL
  )
  
  # variables example: c("count-abandon", "time-of-last-abandon")
  # where.clauses example: c("`product-cart-abandon.item-count` >= 2", "`purchase.total-amount` > 100")
  # extra.events example: c("product-cart-abandon", "page-view", "purchase")
  
  if (length(this$event.attribute)) {
    if (!switch(this$cardinality, 'using.all.values'=, 'using.earliest.value'=, 'using.most.recent.value'=TRUE, FALSE)) {
      stop("Unknown cardinality for event attribute based focalpoint query.  '", this$cardinality, 
           "'.  Must be 'using.all.values', 'using.earliest.value' or 'using.most.recent.value'.")
    }
  } else {
    if (!switch(this$cardinality, 'using.all.events'=, 'using.oldest.event'=, 'using.newest.event'=TRUE, FALSE)) {
      stop("Unknown cardinality for event based focalpoint query.  '", this$cardinality, 
           "'.  Must be 'using.all.events', 'using.oldest.event' or 'using.newest.event'.")
    }
  }
  
  class(this) <- "FocalPointQuery"
  this
}

is.FocalPointQuery <- function(obj) inherits(obj, "FocalPointQuery")

#as.character <- function(this, ...) {
#  UseMethod("as.character", this)
#}

as.character.FocalPointQuery <- function(x, ...) {
  
  # This generates a SQL statement to execute this focalpoint query.
  # The stages are:
  #   Declare the tables we're selecting from, using the focalpoint.event and extra.events. (no dupes)
  #     SELECT * FROM Scenarios variable, Purchase purchase
  #   Declare the start of the where to define the focalpoint point:
  #     WHERE variable.focal_point = purchase.timestamp
  #   Append where.clauses 
  #     AND things from where.clauses
  #   DON'T APPEND LIMIT
  
  # First ensure the focalpoint event is not in the extra events
  #
  events <- unique(x$extra.events)
  if (length(which(events == x$focalpoint.event))) {
    events <- events[-which(events == x$focalpoint.event)]
  }
  x$extra.events <- events
  
  paste(
      "SELECT", GetVariablesSelectClause(x),
    "\n  FROM", GetTablesSelectClause(x), 
    "\n WHERE", GetFocalPointWhereClause(x),
    GetWhereAndClauses(x, "\n   AND"),
    GetLimitClause(x),
    "\n"
  )
}

Limit.FocalPointQuery <- function(this, ...) {
  this$limit
}

`Limit<-.FocalPointQuery` <- function(this, value) {
  this$limit <- value
  this
}

Variables.FocalPointQuery <- function(this, ...) {
  this$variables
}

`Variables<-.FocalPointQuery` <- function(this, value) {
  this$variables <- value
  this
}

# private
# return value examples:
# "`purchase-count`,`page-view-count`"
# "*"
GetVariablesSelectClause <- function(this) {
  stopifnot(is.FocalPointQuery(this))
  if (length(this$variables) == 0) "*" else paste(Backtick(this$variables), collapse=",")
}

# private
# return value examples:
# "Scenarios variables, `focalpoint-event`"
# "Scenarios variables,`extra-event-1`,`cart-abandon`"
GetTablesSelectClause <- function(this) {
  stopifnot(is.FocalPointQuery(this))
  paste(c("Scenarios variable", Backtick(c(this$focalpoint.event, this$extra.events))), collapse=", ")
}

# private
# example return values:
# "variable.focal_point = `purchase`.timestamp"
# "variable.focal_point = `focalpoint-event`.timestamp"
#
GetFocalPointWhereClause <- function(this) {
  stopifnot(is.FocalPointQuery(this))
  paste("variable.focal_point = ", Backtick(this$focalpoint.event), ".", GetFocalPointAttribute(this), sep="")
}

GetFocalPointAttribute <- function(this) {
  if (length(this$event.attribute)) Backtick(this$event.attribute) else "timestamp"
}

# private
# Gets the ANDs that go after the WHERE that specifies the focal point.  These come from this$where.clauses
# 
GetWhereAndClauses <- function(this, separator) {
  sep <- paste(separator, " ", sep="")
  clauses <- c(GetCardinalityClause(this), this$where.clauses)
  if (!length(clauses)) "" else paste(sep, paste(clauses, collapse=sep), sep="")
}

GetCardinalityClause <- function(this) {
  stopifnot(is.FocalPointQuery(this))
  event.bt <- Backtick(this$focalpoint.event)
  
  switch (this$cardinality,
    "using.all.events"=NULL,
    "using.all.values"=NULL,
    "using.oldest.event"=paste("IS_LAST(",  event.bt, ".timestamp)", sep=""),
    "using.newest.event"=paste("IS_FIRST(", event.bt, ".timestamp)", sep=""),
    stop("Unhandled cardinality constraint: ", this$cardinality)
  )
}

# private
# returns NULL or LIMIT XXX
GetLimitClause <- function(this) {
  if (length(this$limit)) paste("LIMIT", this$limit) else NULL
}

AddToFocalPointQuery <- function(focalpoint.query, addition) {
  if (is.WithVariables(addition)) {
    AddWithVariablesToFocalPointQuery(focalpoint.query, addition)
  } else if (is.Where(addition)) {
    AddWhereToFocalPointQuery(focalpoint.query, addition)
  } else if (is.WithEvents(addition)) {
    AddWithEventsToFocalPointQuery(focalpoint.query, addition)
  } else if (is.Limit(addition)) {
    AddLimitToFocalPointQuery(focalpoint.query, addition)
  } else {
    stop("Cannot add ", class(addition), " objects to a FocalPointQuery")
  }
}

AddWithVariablesToFocalPointQuery <- function(this, sql.variables) {
  stopifnot(is.WithVariables(sql.variables))
  this$variables <- sql.variables$variables
  this
}

AddWhereToFocalPointQuery <- function(this, where) {
  stopifnot(is.Where(where))
  this$where.clauses <- c(this$where.clauses, where$clause)
  this
}

AddWithEventsToFocalPointQuery <- function(this, events) {
  stopifnot(is.WithEvents(events))
  this$extra.events <- unique(c(this$extra.events, events$events))
  this
}

AddLimitToFocalPointQuery <- function(this, limit) {
  stopifnot(is.Limit(limit))
  this$limit <- limit$limit
  this
}

WithEvents <- function(...) {
  #
  # TODO: check that these events exist in the config.  (Also ensure the case is correct)
  #
  self <- list(events=as.character(c(list(...), recursive=TRUE)))
  class(self) <- "WithEvents"
  self
}

is.WithEvents <- function(obj) inherits(obj, "WithEvents")


Ops.FocalPointQuery <- function(e1, e2) {
  lhs <- e1 # copy
  rhs <- e2
  if (nargs() == 1) stop("Unary ", .Generic, " not defined for FocalPointQuery")
  if (.Generic != "+") stop("Only + operator is defined for FocalPointQuery")
  
  if (is.FocalPointQuery(lhs)) {
    AddToFocalPointQuery(lhs, rhs)
  } else if (is.WithVariables(rhs)) {
    AddToFocalPointQuery(rhs, lhs)
  } else if (is.Where(rhs)) {
    AddToFocalPointQuery(rhs, lhs)
  } else if (is.WithEvents(rhs)) {
    AddToFocalPointQuery(rhs, lhs)
  } else if (is.Limit(rhs)) {
    AddToFocalPointQuery(rhs, lhs)
  } else {
    NextMethod(.Generic)
  }
}