library("R.methodsS3")

message("TESTING: getDispatchMethodS3()...")

fcn <- getDispatchMethodS3("print", "default")
print(fcn)

tryCatch({
  fcn <- getDispatchMethodS3("print", "unknown")
  print(fcn)
}, error = function(ex) {
  print(ex)
})

message("TESTING: getDispatchMethodS3()...DONE")
