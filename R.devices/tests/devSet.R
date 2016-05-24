message("*** devSet() ...")

library("R.devices")

set.seed(0xBEEF)

devList0 <- devList()

# Open several devices in random order
idxs <- sample(2:20, size=5L, replace=FALSE)
for (idx in idxs) {
  idxT <- devSet(idx)
  # Sanity check
  stopifnot(idxT == idx)
}

# Close all devices
devOff(idxs)

# Sanity check
stopifnot(all.equal(devList(), devList0))


keys <- list(
  A = list(a=1, b=2, c=3),
  B = list(a=1, b=2, c=4),
  C = list(a=1, b=2),
  D = c(a=1, b=2, c=3)
)

# Open several devices using key objects
idxs <- NULL
for (key in keys) {
  idx <- devSet(key)
  idxs <- c(idxs, idx)
}

# Close all devices
devOff(idxs)

# Sanity check
stopifnot(all.equal(devList(), devList0))


# Open several devices using labels
labels <- c("A", "B", "C", "D")

for (rev in c(FALSE, TRUE)) {
  for (label in labels) devSet(label)

  if (rev) labels <- rev(labels)
  for (label in labels) {
    devOff(label)
    left <- devList()
    stopifnot(!is.element(label, names(left)))
  }
}

# Sanity check
stopifnot(all.equal(devList(), devList0))


message("*** devSet() ... DONE")

