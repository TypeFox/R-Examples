## ------------------------------------------------------------------------
library(aoos)

Person <- defineClass("Person", {
  
  personName <- ""
  
  init <- function(name) {
    self$personName(name)
    self$greet()
    }
  
  greet <- function() {
    cat(paste0("Hello, my name is ", self$personName(), ".\n"))
    }
  
  })

## ------------------------------------------------------------------------
ann <- Person("Ann")
ann
ann$personName()
ann$personName("not Ann")
ann$greet()

## ---- eval=FALSE---------------------------------------------------------
#  Person <- defineClass("Person", {
#  
#    .personName <- "" # .personName is private
#  
#    init <- function(name) {
#      self$.personName <- name # option 1
#      .personName <<- name # option 2
#      greet()
#      }
#  
#    greet <- public(function() {
#      cat(paste0("Hello, my name is ", .personName, ".\n")) # before: personName()
#      })
#    })

## ------------------------------------------------------------------------
new("Person", "Ann")

## ------------------------------------------------------------------------
Queue <- defineClass("Queue", {

  .queue <- list()
  
  init <- function(...) {
      for (item in list(...)) self$add(item)
    }
  
  add <- function(x) {
      .queue <<- c(.queue, list(x))
      invisible(self)
    }

  remove <- function() {
      if (queueIsEmpty()) return(NULL)
      head <- .queue[[1]]
      .queue <<- .queue[-1]
      head
    }

  queueIsEmpty <- function() length(.queue) == 0
  
})

HistoryQueue <- defineClass("HistoryQueue", contains = "Queue", {
  
  .head_idx <- 0
  
  remove <- function() {
    if ((length(.queue) - .head_idx) == 0) return(NULL)
    self$.head_idx <- .head_idx + 1
    .queue[[.head_idx]]
    }
  
  show <- function() {
      cat("Next item is at index", .head_idx + 1, "\n")
      for (i in seq_along(.queue)) {
        cat(i, ": ", .queue[[i]], "\n", sep = "")
      }
    }
  
  })

q <- Queue(5, 6, "foo")
q$remove()
q

hq <- HistoryQueue(5, 6, "foo")
hq
hq$show()
hq$remove()
hq$show()
hq$remove()

## ------------------------------------------------------------------------
setMethod("show", "HistoryQueue", function(object) {
  object$show()
})

hq

