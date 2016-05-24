library("R.methodsS3")
library("R.oo")

oopts <- options(warn=1L)

message("TESTING: finalize() on Object on and off ...")

finalized <- NULL
if ("covr" %in% loadedNamespaces()) {
  assertFinalization <- function(name) TRUE
} else {
  assertFinalization <- function(name) {
    cat(sprintf("Is '%s' in '%s'?\n", name, paste(finalized, collapse=", ")))
    stopifnot(is.element(name, finalized))
  }
}

name <- NULL
nextName <- function() {
  if (is.null(name)) return(letters[1L])
  letters[which(letters == name) + 1L]
}

setMethodS3("finalize", "Foo", function(this, ...) {
  cat(sprintf("Finalizing %s()...\n", class(this)[1L]))
  name <- unclass(this)
  cat(sprintf("  Value: %s\n", name))
  finalized <<- c(finalized, name)
  cat(sprintf("Finalizing %s()...done\n", class(this)[1L]))
})


setConstructorS3("Foo", function(..., ...finalize=NA) {
  extend(Object(...), "Foo", ...finalize=...finalize)
})

# Default
x <- Foo(name <- nextName())
rm(list="x"); gc()
assertFinalization(name)

# Default (explicit)
x <- Foo(name <- nextName(), finalize=TRUE, ...finalize=NA)
rm(list="x"); gc()
str(finalized)
assertFinalization(name)

# Disable
x <- Foo(name <- nextName(), finalize=FALSE, ...finalize=FALSE)
rm(list="x"); gc()
str(finalized)

# Disable (forced)
x <- Foo(name <- nextName(), finalize=TRUE, ...finalize=FALSE)
rm(list="x"); gc()
str(finalized)

# Enable (forced)
x <- Foo(name <- nextName(), finalize=FALSE, ...finalize=TRUE)
rm(list="x"); gc()
str(finalized)
assertFinalization(name)

print(finalized)

# Finalize upon exit
options("R.oo::Object/finalizeOnExit"=TRUE)
y <- Foo(name <- "OnExit")

message("TESTING: finalize() on Object on and off ... DONE")

options(oopts)
