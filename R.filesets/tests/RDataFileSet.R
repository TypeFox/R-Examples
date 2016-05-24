library("R.filesets")
loadToEnv <- R.utils::loadToEnv

message("*** RDataFileSet")

x <- list(a=1, b=2)
y <- letters[1:10]
z <- NULL

save(x, file="x.RData")
save(y, file="y.RData")
save(z, file="z.RData")

ds <- RDataFileSet$byPath(".")
print(ds)

x2 <- loadObject(ds[["x"]])
stopifnot(identical(x2, x))

res <- loadObject(ds[["x"]], drop=FALSE)
stopifnot(identical(res$x, x))

y2 <- loadObject(ds[["y"]])
stopifnot(identical(y2, y))

z2 <- loadObject(ds[["z"]])
stopifnot(is.null(z2))

env <- loadToEnv(ds[["x"]])
stopifnot(identical(env$x, x))

message("*** RDataFileSet ... DONE")
