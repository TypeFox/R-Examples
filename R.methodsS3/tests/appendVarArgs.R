library("R.methodsS3")

message("TESTING: appendVarArgs()...")

foobar <- function(a=1) print(a)
print(foobar)

foobar <- appendVarArgs(foobar)
print(foobar)

foobar <- appendVarArgs(foobar)
print(foobar)

# Cleanup
rm(list=ls())

message("TESTING: appendVarArgs()...done")
