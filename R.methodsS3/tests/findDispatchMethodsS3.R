library("R.methodsS3")

message("TESTING: findDispatchMethodS3()...")

## Odds and ends
# Trying to retrieve base::.Options, but should be
# detected as a non-function and return an empty result
fcn <- findDispatchMethodsS3("", "Options")
stopifnot(length(fcn) == 0L)

message("TESTING: findDispatchMethodS3()...DONE")
