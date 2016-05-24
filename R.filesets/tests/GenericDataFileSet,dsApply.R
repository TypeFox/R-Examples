library("R.filesets")
isPackageInstalled <- R.utils::isPackageInstalled
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

message("*** dsApply() on GenericDataFile")

# Example files
path <- system.file("exData", "dataSetA,original", package="R.filesets")
print(path)

# - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up a file set
# - - - - - - - - - - - - - - - - - - - - - - - -
ds <- GenericDataFileSet$byPath(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - -
# Get the size of each file
# - - - - - - - - - - - - - - - - - - - - - - - -
# Alt 1.
message("**** lapply()")
res1 <- lapply(ds, FUN=getFileSize)
str(res1)


# Alt 2. (via an internal loop)
message("**** dsApply(..., .parallel='none')")
res2 <- dsApply(ds, FUN=getFileSize, .parallel="none")
str(res2)
## FIXME: dsApply() returns "short" names whereas lapply() "full" names
stopifnot(all.equal(res2, res1, check.attributes=FALSE))
res1 <- res2 ## FIXME: Workaround trick


# Alt 3a. (via eager futures)
message("**** dsApply(..., .parallel='future') with plan(eager)")
future::plan("eager")
res3a <- dsApply(ds, FUN=getFileSize, .parallel="future")
str(res3a)
stopifnot(all.equal(res3a, res1, check.attributes=FALSE))

# Alt 3b. (via lazy futures)
message("**** dsApply(..., .parallel='future') with plan(lazy)")
future::plan("lazy")
res3b <- dsApply(ds, FUN=getFileSize, .parallel="future")
str(res3b)
stopifnot(all.equal(res3b, res1, check.attributes=FALSE))

# Alt 3c. (via batchjobs futures)
message("**** dsApply(..., .parallel='future') with plan(multicore)")
future::plan("multicore")
res3b <- dsApply(ds, FUN=getFileSize, .parallel="future")
str(res3b)
stopifnot(all.equal(res3b, res1, check.attributes=FALSE))


# Alt 4. (via multicore futures)
if (fullTest && isPackageInstalled("async")) {
  message("**** dsApply(..., .parallel='future') with plan(async::batchjobs, backend='local')")
  ns <- getNamespace("async")
  batchjobs <- get("batchjobs", envir=ns)
  future::plan(batchjobs, backend="local")
  res4 <- dsApply(ds, FUN=getFileSize, .parallel="future")
  str(res4)
  stopifnot(all.equal(res4, res1, check.attributes=FALSE))
}


# Alt 5. (via BatchJobs)
if (fullTest && isPackageInstalled("BatchJobs")) {
  message("**** dsApply(..., .parallel='BatchJobs')")

  ## Override any .BatchJobs.R settings and use a "local" backend
  library("BatchJobs")
  setConfig(list(cluster.functions=makeClusterFunctionsLocal()))

  res4 <- dsApply(ds, FUN=getFileSize, .parallel="BatchJobs")
  print(res4)
  stopifnot(all.equal(res4, res1))
}


# Alt 6. (via BiocParallel + BatchJobs)
if (fullTest && isPackageInstalled("BiocParallel") && isPackageInstalled("BatchJobs")) {
  message("**** dsApply(..., .parallel='BiocParallel::BatchJobs')")

  ## Override any .BatchJobs.R settings and use a "local" backend
  library("BatchJobs")
  setConfig(list(cluster.functions=makeClusterFunctionsLocal()))

  res5 <- dsApply(ds, FUN=getFileSize, .parallel="BiocParallel::BatchJobs")
  print(res5)
  stopifnot(all.equal(res5, res1))
}
