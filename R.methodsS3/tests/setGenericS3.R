library("R.methodsS3")

message("TESTING: setGenericS3()...")

myCat.matrix <- function(..., sep=", ") {
  cat("A matrix:\n")
  cat(..., sep=sep)
  cat("\n")
}

myCat.default <- function(..., sep=", ") {
  cat(..., sep=sep)
  cat("\n")
}

setGenericS3("myCat")

myCat(1:10)
mat <- matrix(1:10, ncol=5)
myCat(mat)

setGenericS3("foo", validators=list(R.methodsS3:::rccValidateSetGenericS3))
setGenericS3("foo<-")

bar.default <- function(...) cat("bar.default\n")
bar <- function(...) cat("bar\n")
setGenericS3("bar")

print(getGenericS3("print"))


# Your defintion will redefine bar() above to bar.default().
foobar <- function() print("foobar()")
setGenericS3("foobar")


# Cleanup
rm(list=ls())

message("TESTING: setGenericS3()...DONE")
