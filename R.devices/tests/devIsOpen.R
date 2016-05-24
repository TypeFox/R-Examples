message("*** devIsOpen() ...")

library("R.devices")

for (which in list(c(), 1L, 1:5)) {
  print(which)
  res <- devIsOpen(which)
  print(res)
  stopifnot(is.logical(res))
  stopifnot(is.character(names(res)))
  stopifnot(length(res) == length(which))
}

message("*** devIsOpen() ... DONE")
