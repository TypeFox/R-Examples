message("*** devOptions() ...")

# Without attaching package
opts0 <- R.devices::devOptions()
print(opts0)

opts0eps <- R.devices::devOptions("eps")
str(opts0eps)

opts <- R.devices::devOptions("png")
print(opts)

# With attaching package
library("R.devices")
opts1 <- R.devices::devOptions()
print(opts1)

opts1eps <- R.devices::devOptions("eps")
str(opts1eps)

stopifnot(identical(opts1eps, opts0eps))
stopifnot(identical(opts1, opts0))


# Options for the PNG device
opts <- devOptions("png")
print(opts)

# Options for the postscript device
opts <- devOptions("postscript")
print(opts)

# Same using alias
opts2 <- devOptions("ps")
print(opts2)
stopifnot(identical(opts2, opts))

# Options for all known devices
opts <- devOptions()
print(opts)

# Setting a custom option
devOptions("png", foo=list(a=1, b=pi))
str(devOptions("png")$foo)

# Setting option to NULL, i.e. drop it
devOptions("png", foo=NULL)
str(devOptions("png")$foo)
str(devOptions("png"))

# Get individual device options
print(getDevOption("png", "width"))

opts1 <- R.devices::devOptions()
print(opts1)

## Assert devOptions(<function>) equals devOptions(<name>)
## Identify all "primitive" functions
fcns <- lapply(R.devices:::devAll(), FUN=`[`, 1L)
for (name in names(fcns)) {
  fcn <- fcns[[name]]
  cat(sprintf("Asserting devOptions('%s') == devOptions(%s) ...\n", name, fcn))
  fcn <- eval(parse(text=fcn))
  optsName <- devOptions(name)
  optsFcn <- devOptions(fcn)
  res <- all.equal(optsFcn, optsName)
  if (!isTRUE(res)) {
    str(list(name=name, fcn=fcn, byName=optsName, byFcn=optsFcn))
    print(res)
    stop(!all.equal(optsFcn, optsName))
  }
}

message("*** devOptions() for each type ...")
types <- rownames(devOptions())
cat("All known types:\n")
print(types)

for (type in types) {
  message(sprintf("*** devOptions('%s') for each type ...", type))
  opts <- devOptions(type)
  str(opts)
  message(sprintf("*** devOptions('%s') for each type ... DONE", type))
}

message("*** devOptions() for each type ... DONE")


message("*** devOptions(drop=FALSE) ...")

opts <- devOptions("png", drop=TRUE)
str(opts)
stopifnot(is.list(opts))
stopifnot(is.null(dim(opts)))

opts <- devOptions("png", drop=FALSE)
str(opts)
stopifnot(is.list(opts))
stopifnot(!is.null(dim(opts)))

message("*** devOptions(drop=FALSE) ... DONE")

message("*** devOptions(reset=TRUE) ...")

## Reset all
opts <- devOptions()
print(opts)
opts0 <- devOptions(reset=TRUE)
print(opts0)

## Reset one device
opts <- devOptions("png")
width <- getDevOption("png", "width")
devOptions("png", width=2*width)
stopifnot(getDevOption("png", "width") == 2*width)
devOptions("png", reset=TRUE)
stopifnot(getDevOption("png", "width") == width)

## Reset "*"
opts <- devOptions("*")
path <- getDevOption("*", "path")
devOptions("*", path="foo")
stopifnot(getDevOption("*", "path") == "foo")
devOptions("*", reset=TRUE)
stopifnot(getDevOption("*", "path") == path)


message("*** devOptions(reset=TRUE) ... DONE")


message("*** devOptions() - errors ...")

res <- try(devOptions(type=character(0L), width=32L))
stopifnot(inherits(res, "try-error"))

message("*** devOptions() - errors ... DONE")


message("*** devOptions() - odds'n'ends ...")

devOptions(type=character(0L), reset=FALSE)
devOptions(type=character(0L), reset=TRUE)

message("*** devOptions() - odds'n'ends ... DONE")

message("*** devOptions() ... DONE")
