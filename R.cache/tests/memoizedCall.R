library("R.cache")

# Use an empty temporary file cache
setCacheRootPath(path=file.path(tempdir()))
clearCache(recursive=TRUE, prompt=FALSE)
dirs <- c("tests", "memoizedCall")

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
# There will be no cache hit for the first call
t0 <- system.time({
  res0 <- memoizedCall(sleep, time=1.5, dirs=dirs)
})[3]
print(t0)

# The second will have a cache hit and therefore
# return the memoized results momentarily.
t1 <- system.time({
  res1 <- memoizedCall(sleep, time=1.5, dirs=dirs)
})[3]
print(t1)
if (t1 >= t0) {
  warning("Second call to memoizedCall() took longer than the first: ",
          t1, " >= ", t0)
}

# Sanity check
stopifnot(identical(res1, res0))
clearCache(recursive=TRUE)
