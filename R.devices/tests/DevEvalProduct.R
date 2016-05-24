view <- R.devices:::view

message("*** DevEvalProduct ...")

message("*** DevEvalProduct - subsetting ...")

p <- R.devices:::DevEvalProduct("foo", tags=c("a", "b"))
for (ff in c("fullname", "name", "tags")) {
  cat(sprintf("%s: %s\n", ff, substring(p[[ff]], 1, 50)))
}

library("R.devices")

p <- DevEvalProduct("foo", tags=c("a", "b"))
for (ff in c("fullname", "name", "tags")) {
  cat(sprintf("%s: %s\n", ff, substring(p[[ff]], 1, 50)))
}

p <- DevEvalProduct("foo", tags=c("a", "b"))
valueA <- p[["name"]]
valueB <- p$name
stopifnot(identical(valueB, valueA))
value <- p[["non-existing-field"]]
stopifnot(is.null(value))

message("*** DevEvalProduct - subsetting ... DONE")


message("*** DevEvalProduct - view() ...")
view(p)
!p
message("*** DevEvalProduct - view() ... DONE")

message("*** DevEvalProduct ... DONE")
