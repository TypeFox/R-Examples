message("TESTING: InternalErrorException...")

library("R.oo")

ex <- InternalErrorException("Hmm... I didn't expect this!")
print(ex)

ex2 <- InternalErrorException("Hmm... I didn't expect this!", package=R.oo)
print(ex2)

ex3 <- InternalErrorException("Hmm... I didn't expect this!", package="R.oo")
print(ex3)


myLog <- function(x, ...) {
  if (!is.numeric(x)) {
    throw(InternalErrorException("Argument 'x' to myLog() is not numeric: ",
                                                 mode(x), package=R.oo))
  }
  log(x, ...)
}


myLog(2)

ex <- NULL
tryCatch({
  myLog("a")
}, error= function(ex) {
  ex <- Exception$getLastException()
})

message("TESTING: InternalErrorException...DONE")
