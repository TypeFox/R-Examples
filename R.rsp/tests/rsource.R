library("R.rsp")

rm(list="i")

pathname <- system.file("exData", "slowcounting.txt.rsp", package="R.rsp")
rsource(pathname)

# Assert that the evaluation was done in the current environment
stopifnot(exists("i", mode="numeric"))
