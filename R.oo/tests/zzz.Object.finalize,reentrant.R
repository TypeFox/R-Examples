message("TESTING: finalize() reentrant...")

## FIXME: 'covr' does not play well with tests
## detaching/unloading packages
if ("covr" %in% loadedNamespaces()) {
  detach <- function(...) NULL
}


library("R.methodsS3")
library("R.oo")

lotsOfParsing <- function(code="y <- 1:3") {
  parse(text=rep(code, times=10000))
}

setConstructorS3("MyClass", function(a=1:10) {
  extend(Object(), "MyClass", a=a)
})

setMethodS3("finalize", "MyClass", function(this, ...) {
  cat("finalize...\n")
  utils::str(sys.calls())
  cat("finalize...done\n")
})

# Parse and eval expression (works)
expr <- lotsOfParsing()
eval(expr)
print(y)
## [1] 1 2 3
stopifnot(identical(y, 1:3))

# Create an object with a finalizer
x <- MyClass()

# Detach R.oo so that the finalizer will try to reload it
detach("package:R.oo")

# Remove 'x' so that it will be finalized below
rm(x)


# This may trigger garbage collection via parse()
# If so, it is important that parse() is not called
# (indirectly via library()) by the finalizer.
# Because otherwise R may crash.
expr2 <- lotsOfParsing(code="y <- 1:4")
## finalize...
## Dotted pair list of 9
##  $ : ...
##  ...
##  $ : language function (env)  { ...
##  $ : language finalize(this)
##  $ : language finalize.MyClass(this)
## Parse called: TRUE
## finalize...done
eval(expr2)
print(y)
## [1] 1 2 3 4
stopifnot(identical(y, 1:4))

print(warnings())

message("TESTING: finalize() reentrant...DONE")
