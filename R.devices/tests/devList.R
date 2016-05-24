message("*** devList() ...")

library("R.devices")

res <- devList()
print(res)
stopifnot(is.integer(res))
stopifnot(is.character(names(res)))

res <- devList(dropNull=FALSE)
print(res)
stopifnot(is.integer(res))
stopifnot(is.character(names(res)))

res <- devList(interactiveOnly=TRUE)
print(res)
stopifnot(is.integer(res))
stopifnot(is.character(names(res)))

# - - - - - - - - - - - - - - - - - - - - - - - - -
# Labels
# - - - - - - - - - - - - - - - - - - - - - - - - -
devSetLabel(which=1L, label="foo")
label <- devGetLabel(1L)
print(label)
stopifnot(label == "foo")
label <- devGetLabel("foo")
print(label)
stopifnot(label == "foo")

devSetLabel(which=1L, label="bar")
label <- devGetLabel(1L)
print(label)
stopifnot(label == "bar")
label <- devGetLabel("bar")
print(label)
stopifnot(label == "bar")

devSetLabel(which="bar", label="foo")
label <- devGetLabel(1L)
print(label)
stopifnot(label == "foo")
label <- devGetLabel("foo")
print(label)
stopifnot(label == "foo")


res <- try(devGetLabel(which=10L))
stopifnot(inherits(res, "try-error"))

message("*** devList() ... DONE")

