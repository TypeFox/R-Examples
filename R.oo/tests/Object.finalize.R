library("R.oo")

oopts <- options(warn=1L)

message("TESTING: finalize() on Object ...")

setConstructorS3("MyClass", function() {
  extend(Object(), "MyClass")
})

setMethodS3("finalize", "MyClass", function(this, ...) {
  cat(as.character(this), "is about to be removed from the memory!\n")
})

o <- MyClass()
o <- MyClass()
o <- MyClass()
o <- MyClass()
gc()

## MyClass: 0x01BE602C is about to be removed from the memory!
## MyClass: 0x01BFF634 is about to be removed from the memory!
## MyClass: 0x01C13584 is about to be removed from the memory!
##          used (Mb) gc trigger (Mb)
## Ncells 229903  6.2     467875 12.5
## Vcells  53725  0.5     786432  6.0

rm(o)
## MyClass: 0x01C578B0 is about to be removed from the memory!
##          used (Mb) gc trigger (Mb)
## Ncells 229903  6.1     467875 12.3
## Vcells  53725  0.5     786432  6.0

message("TESTING: finalize() on Object ... DONE")

options(oopts)
