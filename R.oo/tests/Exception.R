message("TESTING: Exception...")

library("R.oo")
oopts <- options(warn=1)

######################################################################
# 1. To catch a regular "error" exception thrown by e.g. stop().
######################################################################
x <- NA
y <- NA
tryCatch({
  x <- log(123)
  y <- log("a")
}, error = function(ex) {
  print(ex)
})
print(x)
print(y)



######################################################################
# 2. Always run a "final" expression regardless or error or not.
######################################################################
filename <- tempfile("R.methodsS3.example")
con <- file(filename)
tryCatch({
  open(con, "r")
}, error = function(ex) {
  cat("Could not open ", filename, " for reading.\n", sep="")
}, finally = {
  close(con)
  cat("The id of the connection is ",
       ifelse(is.null(con), "NULL", con), ".\n", sep="")
})


######################################################################
# 3. Creating your own Exception class
######################################################################
setConstructorS3("NegativeLogValueException", function(
  msg="Trying to calculate the logarithm of a negative value", value=NULL) {
  extend(Exception(msg=msg), "NegativeLogValueException",
    .value = value
  )
})

setMethodS3("as.character", "NegativeLogValueException", function(this, ...) {
  paste(as.character.Exception(this), ": ", getValue(this), sep="")
})

setMethodS3("getValue", "NegativeLogValueException", function(this, ...) {
  this$.value
})


mylog <- function(x, base=exp(1)) {
  if (x < 0)
    throw(NegativeLogValueException(value=x))
  else
    log(x, base=base)
}


# Note that the order of the catch list is important:
l <- NA
x <- 123
tryCatch({
  l <- mylog(x)
}, NegativeLogValueException = function(ex) {
  cat(as.character(ex), "\n")
}, "try-error" = function(ex) {
  cat("try-error: Could not calculate the logarithm of ", x, ".\n", sep="")
}, error = function(ex) {
  cat("error: Could not calculate the logarithm of ", x, ".\n", sep="")
})
cat("The logarithm of ", x, " is ", l, ".\n\n", sep="")

options(oopts)

message("TESTING: Exception...DONE")
