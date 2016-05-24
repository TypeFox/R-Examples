library("R.cache")

# Use an empty temporary file cache
setCacheRootPath(path=file.path(tempdir()))
clearCache(recursive=TRUE, prompt=FALSE)
dirs <- c("tests", "addMemoization")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define function to be memoized
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sleep <- function(time) {
  cat(sprintf("Sleeping for %g seconds...\n", time))
  Sys.sleep(time)
  cat(sprintf("Sleeping for %g seconds...done\n", time))
  time
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test memoization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sleep <- addMemoization(sleep)

# There will be no cache hit for the first call
t0 <- system.time({
  res0 <- sleep(1.5)
})[3]
print(t0)

# The second will have a cache hit and therefore
# return the memoized results momentarily.
t1 <- system.time({
  res1 <- sleep(1.5)
})[3]
print(t1)
print(t1/t0)

# Sanity check
stopifnot(identical(res1, res0))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Don't memoize already memoized functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sleep2 <- addMemoization(sleep)
stopifnot(identical(sleep2, sleep))

sleep3 <- addMemoization("sleep")
stopifnot(identical(sleep3, sleep))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exception handling
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- try(addMemoization("non-existing-function"), silent=TRUE)
stopifnot(inherits(res, "try-error"))
