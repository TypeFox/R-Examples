library("R.cache")
oopts <- options("R.cache::compress")

simulate <- function(mean, sd) {
  # 1. Try to load cached data, if already generated
  key <- list(mean, sd)
  data <- loadCache(key)
  if (!is.null(data)) {
    cat("Loaded cached data\n")
    return(data)
  }

  # 2. If not available, generate it.
  cat("Generating data from scratch...")
  data <- rnorm(1000, mean=mean, sd=sd)
  Sys.sleep(1)             # Emulate slow algorithm
  cat("ok\n")
  saveCache(data, key=key, comment="simulate()")

  data
}


for (compress in c(FALSE, TRUE)) {
  options("R.cache::compress"=compress)

  data <- simulate(2.3, 3.0)
  data <- simulate(2.3, 3.5)
  data <- simulate(2.3, 3.0) # Will load cached data

  # Clean up
  file.remove(findCache(key=list(2.3,3.0)))
  file.remove(findCache(key=list(2.3,3.5)))
}

## Cleanup
options(oopts)
