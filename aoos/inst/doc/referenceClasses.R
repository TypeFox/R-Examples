## ------------------------------------------------------------------------
library(aoos)

Person <- defineRefClass({
  Class <- "person" # this is the argument 'Class' in setRefClass
  
  personName <- "character" # this is a field of class 'character'
  
  initialize <- function(name) {
    .self$personName <- name
    .self$greet()
  }
  
  greet <- function() {
    cat(paste0("Hello, my name is ", .self$personName, ".\n"))
  }
  
})

ann <- Person("Ann")
ann
ann$personName
ann$personName <- "not Ann"
ann$greet()

## ------------------------------------------------------------------------
PrivatePerson <- defineRefClass({
  Class <- "PrivatePerson"
  contains <- "Private" # also just passed as argument to setRefClass
  
  .personName <- "character"
  
  initialize <- function(name) {
    .self$.personName <- name
    .self$greet()
  }
  
  greet <- function() {
    cat(paste0("Hello, my name is ", .self$.personName, ".\n"))
  }
  
})

ann <- PrivatePerson("Ann")
ann
stopifnot(inherits(try(ann$.personName, silent = TRUE), "try-error"))
ann$greet()

## ------------------------------------------------------------------------
removeClass("PrivatePerson")

PrivatePerson <- setRefClass(
  Class = "PrivatePerson", 
  fields = list(.personName = "character"),
  contains = "Private",
  methods = list(
    initialize = function(name) {
      .self$.personName <- name
      .self$greet()
    },
    greet = function() {
      cat(paste0("Hello, my name is ", .self$.personName, ".\n"))
    }
  )
)

ann <- PrivatePerson("Ann")
ann
stopifnot(inherits(try(ann$.personName, silent = TRUE), "try-error"))
ann$greet()

