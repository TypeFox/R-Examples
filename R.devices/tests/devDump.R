message("*** devDump() ...")

library("R.devices")
graphics.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Copy content of current screen device
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  # Open device
  plot(1:10)

  devList0 <- devList()
  devEval("png,jpg,pdf", name="count", tags="copy")
  stopifnot(all.equal(devList(), devList0))

  # Sanity checks
  print(devList())
  stopifnot(length(devList()) == 1L)


  # Same using a default name
  devList0 <- devList()
  devEval("png,jpg,pdf")
  stopifnot(all.equal(devList(), devList0))
  stopifnot(length(devList()) == 1L)

  # Close device
  devOff()

  stopifnot(length(devList()) == 0L)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy content of all screen devices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open several devices
  idxs <- NULL

  idxs <- c(idxs, devNew())
  plot(1:10)

  idxs <- c(idxs, devNew())
  plot(cos)

  # Automatially "dump" image files of all open devices
  devDump()

  # "Manual" saving of defined
  which <- devList()
  print(which)

  # Save all
  devEval("png,pdf", which=which)

  # Close all opened devices
  devOff(idxs)
} # if (interactive())


message("*** devDump() ... DONE")
