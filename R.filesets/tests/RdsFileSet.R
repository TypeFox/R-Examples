library("R.filesets")

message("*** RdsFileSet")

x <- list(a=1, b=2)
y <- letters[1:10]

saveRDS(x, "x.rds")
saveRDS(y, "y.rds")

ds <- RdsFileSet$byPath(".")
print(ds)

x2 <- loadObject(ds[[1]])
stopifnot(identical(x2, x))

message("*** RdsFileSet ... DONE")
