message("TESTING: objectSize()...")

library("R.oo")

## Simple object
x <- 1:100
y0 <- object.size(x)
y <- objectSize(x)
print(y)
stopifnot(y == y0)

## A list
x <- as.list(1:100)
y0 <- object.size(x)
y <- objectSize(x)
print(y)
stopifnot(y == y0)

## An environment
env <- new.env()
env$a <- 1:100
env$b <- as.list(1:100)
env$c <- new.env()
y0 <- object.size(env)
print(y0)
y <- objectSize(env)
print(y)

## An environment with circular dependencies
env <- new.env()
env$a <- 1:100
env$env <- env
y0 <- object.size(env)
print(y0)
y <- objectSize(env)
print(y)

message("TESTING: objectSize()...DONE")
