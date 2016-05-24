message("TESTING: throw()...")

library("R.oo")

## Generate an error
ex <- tryCatch({
  stop("An error")
}, error = function(ex) {
  ex
})
print(ex)

## Re-throw the error
ex2 <- tryCatch({
  throw(ex)
}, error = function(ex) {
  ex
})
print(ex2)

stopifnot(identical(ex2, ex))


## Generate an Exception
ex <- tryCatch({
  throw("An error")
}, error = function(ex) {
  ex
})
print(ex)

## Re-throw the Exception
ex2 <- tryCatch({
  throw(ex)
}, error = function(ex) {
  ex
})
print(ex2)

stopifnot(identical(ex2, ex))


message("TESTING: throw()...DONE")
