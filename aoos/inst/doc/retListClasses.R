## ----eval=FALSE----------------------------------------------------------
#  constructor <- function(...) {
#    ...
#    out <- list(...)
#    class(out) <- "constructor"
#    out
#  }

## ----eval=FALSE----------------------------------------------------------
#  constructor <- function(...) {
#    ...
#    retList("constructor")
#  }

## ------------------------------------------------------------------------
library("aoos")

Person <- function(name) {
  
  print <- function(x, ...) {
    cat(paste0("Hello, my name is ", .self$name, ".\n"))
  }
  
  retList(c("Person", "Print"))
}

## ------------------------------------------------------------------------
ann <- Person("Ann")
ann

## ------------------------------------------------------------------------
Person <- function(.name) {
  
  print <- function(x, ...) {
    cat(paste0("Hello, my name is ", .self$.name, ".\n"))
  }
  
  name <- function(x) {
    if (!missing(x)) .name <<- x
    .name
  }
  
  retList(c("Person", "Print"))
}

p <- Person("Ann")
p
p$name()
p$name("Paul")
p

## ------------------------------------------------------------------------
Person <- function(name) {
  
  print <- function(x, ...) {
    cat(paste0("Hello, my name is ", .self$name, ".\n"))
  }
  
  retList(c("Person", "Print"))

}

Employee <- function(id, ...) {
  
  print <- function(x, ...) {
    cat("Name: ", name, "\nID:   ", id)
  }
  
  retList("Employee", super = Person(...))
  
}

kalle <- Employee("1", "Kalle")
str(kalle)
kalle

