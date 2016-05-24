message("TESTING: getConstructorS3()...")

library("R.oo")


fcn <- getConstructorS3("Object")
print(args(fcn))

message("TESTING: getConstructorS3()...DONE")
