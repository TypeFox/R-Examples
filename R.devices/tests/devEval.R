message("*** devEval() ...")

library("R.devices")
library("R.utils")
graphics.off()

png <- grDevices::png
jpeg <- grDevices::jpeg
tiff <- grDevices::tiff


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Various types of single and multiple device outputs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval() - single and multiple device outputs ...")

types <- list(
  character(0L),
  "png",
  "jpg",
  c("png", "png", "jpeg"),
  "png,jpg,pdf"
)

for (type in types) {
  cat("Device types: ", paste(sQuote(type), collapse=", "), "\n", sep="")
  devList0 <- devList()
  res <- devEval(type, name="multi", aspectRatio=2/3, {
    plot(1:10)
  })
  print(res)
  stopifnot(length(res) == length(unlist(strsplit(type, split=","))))
  stopifnot(all.equal(devList(), devList0))
}

# Sanity checks
print(devList())
stopifnot(length(devList()) == 0L)

message("*** devEval() - single and multiple device outputs ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# With 'initially' and 'finally' expression
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval() - initially and finally ...")

devList0 <- devList()
devEval(c("png", "jpg"), name="count", {
  plot(1:10)
  count <- count + 1L
}, initially = {
  # Emulate an overhead
  cat("Initiate...")
  count <- 0L
  Sys.sleep(1)
  cat("done\n")
}, finally = {
  cat("Number of image files created: ", count, "\n", sep="")
})
stopifnot(all.equal(devList(), devList0))

# Sanity checks
print(devList())
stopifnot(length(devList()) == 0L)

message("*** devEval() - initially and finally ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Try several devices until first successful device is found
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval() - first successful device ...")

types <- list(
  "png|jpg|pdf",               # PNG, JPG and then PDF
  "dummy|png|jpg|pdf",         # "Non-existing", PNG, JPG and then PDF
  "quartz|x11|windows",        # Any interactive device (depending on OS)
  c("png|jpg", "x11|windows"), # PNG or JPG and then x11 or windows
  "eps|postscript|pdf",        # EPS, Postscript or PDF
  "jpeg2|jpeg",                # JPEG via bitmap() or via jpeg()
  "png,jpg|x11|windows"        # == c("png", "jpg|x11|windows")
)

if (!capabilitiesX11()) {
  message("Skipping test for X11")
  types <- lapply(types, FUN=function(x) gsub("x11|", "", x, fixed=TRUE))
}

devList0 <- devList()

for (type in types) {
  printf("Any of %s\n", paste(sQuote(type), collapse=" + "))

  # Use try-catch in case not supported on some test systems
  tryCatch({
    res <- devEval(type, name="any", aspectRatio=2/3, scale=1.2, {
      plot(100:1)
    })
    printf("Result: %s (%s)\n\n", sQuote(res), attr(res, "type"))

    if (length(devList()) > 0) devOff()
  }, error = function(ex) {
    printf("Failed: %s\n\n", sQuote(ex$message))
  })
} # for (type ...)

# Sanity check
stopifnot(all.equal(devList(), devList0))

message("*** devEval() - first successful device ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot a parsed expression
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval() - parsed expressions ...")

expr <- substitute(plot(1:10))
tryCatch({
  res <- devEval("png|jpg", name="any", width=480L, height=480L, {
    plot(100:1)
  })
  printf("Result: %s (%s)\n\n", sQuote(res), attr(res, "type"))

  if (length(devList()) > 0) devOff()
}, error = function(ex) {
  printf("Failed: %s\n\n", sQuote(ex$message))
})

message("*** devEval() - parsed expressions ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Special cases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** toDefault(<expr>) ...")

# toX11({ plot(1:10) }) actually results in a call to
# devEval(type="x11", name={ plot(1:10) }); note argument 'name'
# and not 'expr'.  The following tests that devEval() recognizes
# and handles this internally.

## FIXME: The current solution evaluates 'name' internally
## and therefore opens a interactive graphics device.
if (interactive()) {
   res <- toDefault({ plot(1:10) })
   print(res)

   ## FIX ME:
   graphics.off()
}

message("*** toDefault(<expr>) ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Device type specified as a device functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval(<fcn>) ...")

types <- list(
  png=grDevices::png,
  jpg=grDevices::jpeg
)

for (name in names(types)) {
  cat("Device types: ", paste(sQuote(name), collapse=", "), "\n", sep="")
  type <- types[[name]]
  str(args(type))
  devList0 <- devList()
  res <- devEval(type, name="multi", tags="function", aspectRatio=2/3, {
    plot(1:10)
  })
  print(res)
  stopifnot(length(res) == length(type))
  stopifnot(all.equal(devList(), devList0))
}

message("*** devEval(<fcn>) ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Special case: Default device
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** devEval(<default>) ...")

cat("Device types: 'default'\n")
type <- getOption("device")
str(type)
devList0 <- devList()
res <- devEval(type, name="default", aspectRatio=2/3, {
  plot(1:10)
})
print(res)
wasInteractiveOpened <- (length(setdiff(devList(), devList0)) > 0L)
if (wasInteractiveOpened) devOff()

message("*** devEval(<default>) ... DONE")



# Sanity checks
print(devList())
stopifnot(length(devList()) == 0L)

message("*** devEval() ... DONE")
