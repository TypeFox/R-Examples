## ------------------------------------------------------------------------
library("R6")
library("aoos")
library("microbenchmark")

R6 <- R6Class("R6",
              public = list(
                x = NULL,
                initialize = function(x = 1) self$x <- x,
                getx = function() self$x,
                inc = function(n = 1) self$x <- x + n
              )
)

RC <- setRefClass("RC", 
                  fields = list(x = "numeric"),
                  methods = list(
                    initialize = function(x = 1) .self$x <- x,
                    getx = function() x,
                    inc = function(n = 1) x <<- x + n
                  )
)

RList <- function(x = 1) {
  self <- environment()
  getx <- function() self$x
  inc <- function(n = 1) self$x <- self$x + n
  out <- list(x = x, getx = getx, inc = inc)
  class(out) <- "RList"
  out
}

RL <- function(x = 1) {
  getx <- function() .self$x
  inc <- function(n = 1) .self$x <- .self$x + n
  retList("RL", c("x", "getx", "inc"))
}

DC <- defineClass("DC", {
  getx <- function() .self$x
  inc <- function(n = 1) .self$x <- .self$x + n
  init <- function(x = 1) .self$x <- x
})

# And some more definitions for inheritance
R6Child <- R6Class("R6Child", inherit = R6)
RCChild <- setRefClass("RCChild", contains = "RC")
RLChild <- function(...) {
  retList("RLChild", super = RL(...))
}
DCChild <- defineClass("DCChild", contains = "DC", {})

## ------------------------------------------------------------------------
benchmark1 <- microbenchmark(
  DC(),
  RC$new(),
  R6$new(),
  RL(),
  RList()
) 

print(benchmark1, digits = 2)

## ------------------------------------------------------------------------
benchmark2 <- microbenchmark(
  DCChild(),
  RCChild$new(),
  R6Child$new(),
  RLChild()
)

print(benchmark2, digits = 2)

