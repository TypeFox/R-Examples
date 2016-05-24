library("R.methodsS3")

message("TESTING: setMethodS3()...")

######################################################################
# Example 1
######################################################################
setMethodS3("foo", "default", function(x, ...) {
  cat("In default foo():\n")
  print(x, ...)
})


setMethodS3("foo", "character", function(s) {
  cat("In foo() for class 'character':\n")
  print(s, ...)
})

# The generic function is automatically created!
print(foo)

foo(123)
foo("123")


######################################################################
# Example 2
#
# Assume that in a loaded package there is already a function bar(),
# but you also want to use the name 'bar' for the character string.
# It may even be the case that you do not know of the other package,
# but your users do!
######################################################################
# bar() in other package
bar <- function(x, y, ...) {
  cat("In bar() of 'other' package.\n")
}


# Your defintion will redefine bar() above to bar.default().
setMethodS3("bar", "character", function(object, ...) {
  cat("In bar() for class 'character':\n")
  print(object, ...)
})

bar(123)
bar("123")

setMethodS3("bar<-", "character", function(x, value) {
  attr(x, "bar") <- value
  x
})

x <- "a"
bar(x) <- "hello"
str(x)


setMethodS3("$", "SomeClass", function(x, name) {
  attr(x, name)
})

setMethodS3("$<-", "SomeClass", function(x, name, value) {
  attr(x, name) <- value
  x
})



setMethodS3("yaa", "character", abstract=TRUE, validators=list(R.methodsS3:::rccValidateSetMethodS3))

print(getMethodS3("yaa", "character"))

# Redefine
setMethodS3("yaa", "character", abstract=TRUE, validators=list(R.methodsS3:::rccValidateSetMethodS3))



# Cleanup
rm(list=ls())

message("TESTING: setMethodS3()...DONE")
