message("*** devEval() - errors ...")

library("R.devices")
hpaste <- R.utils::hpaste
graphics.off()

path <- getDevOption("png", "path")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Error handling: Incomplete image file is removed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- getDevOption("png", "path")
files0 <- dir(path=path)

tryCatch({
  res <- devEval(type="png", name="error", {
    plot(1:10)
    v <- log("0")
    abline(v=v)
  }, onIncomplete="remove")
  print(res)
}, error = function(ex) {
  message("An error occurred while plotting: ", ex$message)
})

## Assert that any image created was removed
files <- dir(path=path)
new <- setdiff(files, files0)
if (length(new) > 0L) {
  stop("Failed to remove incomplete image file: ", hpaste(sQuote(new)))
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Error handling: Incomplete image file is renamed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- getDevOption("png", "path")

## Try many times to test unique renaming of files.
for (kk in 1:5) {
  files0 <- dir(path=path)
  tryCatch({
    res <- devEval(type="png", name="error", {
      plot(1:10)
      v <- log("0")
      abline(v=v)
    }, onIncomplete="rename")
    print(res)
  }, error = function(ex) {
    message("An error occurred while plotting: ", ex$message)
  })

  ## Assert that any image created was removed
  files <- dir(path=path)
  new <- setdiff(files, files0)
  if (length(new) != 1L) {
    stop("Failed to rename incomplete image file")
  }
  cat("Incomplete image file: ", new, "\n", sep="")
}

# Sanity checks
print(devList())
stopifnot(length(devList()) == 0L)

message("*** devEval() - errors ... DONE")
