message("*** devListIndexOf() ...")

library("R.devices")

.devListIndexOf <- R.devices:::.devListIndexOf

labels <- list(character(0L), "Device 1", c("Device 1", "Device 1"), c("Device 1", "FooBar"))
for (label in labels) {
  res <- .devListIndexOf(label, error=FALSE)
  print(res)
  stopifnot(is.integer(res));
  stopifnot(is.character(names(res)));
}

message("*** devListIndexOf() ... DONE")
