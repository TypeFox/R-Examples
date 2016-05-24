oopts <- options(warn=1L)

message("TESTING: finalize() on Object without attach ...")

pkgs <- c("R.methodsS3", "R.oo")
isAttached <- function(pkgs) {
  structure(sprintf("package:%s", pkgs) %in% search(), names=pkgs)
}

# Record which packages were attached from the beginning
# (happens if not a fresh session)
wasAttached <- isAttached(pkgs)


assertPackages <- function(loaded=c("R.methodsS3", "R.oo")) {
  s <- utils::sessionInfo()
  s$R.version <- NULL;
  s$platform <- "";
  s$locale <- "";
  cat("----------------------------------")
  print(s)
  cat("----------------------------------\n\n")
  loaded <- loaded[!wasAttached[loaded]]
  stopifnot(!any(isAttached(loaded)))
}

R.oo::setConstructorS3("MyClass", function(a=1:10) {
  R.oo::extend(R.oo::Object(), "MyClass", a=a)
})

# Create an object with a finalizer
x <- MyClass()

assertPackages()

# Remove 'x' so that it will be finalized below
rm(x)
gc()

assertPackages(loaded="R.oo")

message("TESTING: finalize() on Object without attach ... DONE")

options(oopts)
