message("TESTING: ll()...")

library("R.oo")

## Create some objects in the current environment
a <- 1:100
env <- new.env()
env$b <- letters[1:10]

## List the content of evironments
x <- ll(envir=env)
print(x)

## Empty environment
x <- ll(envir=new.env())
print(x)

## search() path environment
x <- ll(envir=1L)
str(x)

## search() path environment
x <- ll(envir="R.oo")
str(x)

## Filter by name pattern
x <- ll(envir="R.oo", pattern="^throw.*")
print(x)

x <- ll(envir="R.oo", pattern="^NonExistingName$")
print(x)

## List all functions and sort them by size
x <- ll(envir="R.oo", mode="function", sortBy="objectSize")
str(x)

## List all functions of a package and sort them by size
x <- ll(R.oo, mode="function", sortBy="objectSize")
str(x)

## List all functions of an Object
obj <- Object()
obj$a <- 1:100
obj$b <- new.env()
x <- ll(obj)
print(x)

message("TESTING: ll()...DONE")
