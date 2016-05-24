message("TESTING: attachLocally()...")

library("R.oo")

foo <- function(object, arg1="some value", ...) {
  cat("Local objects in foo():\n")
  print(ls())

  attachLocally(object)

  cat("\nLocal objects in foo():\n")
  print(ls())

  for (name in ls()) {
    cat("\nObject '", name, "':\n", sep="")
    print(get(name, inherits=FALSE))
  }
}

a <- "A string"
obj <- Object()
obj$a <- "Another string"
obj$b <- NA
foo(obj)
print(a)

message("TESTING: attachLocally()...DONE")
