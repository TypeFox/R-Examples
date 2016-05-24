library("R.devices")

message("*** toNnn() ...")

## Get all toNnn() functions
envir <- getNamespace("R.devices")
names <- ls(pattern="^to", envir=envir)
fcns <- list()
for (name in names) {
  if (exists(name, envir=envir, mode="function")) {
    fcns[[name]] <- get(name, envir=envir, mode="function")
  }
}

for (name in names(fcns)) {
  message(sprintf("*** %s() ...", name))
  toFcn <- fcns[[name]]
  tryCatch({
    toFcn(name=name, { plot(1:10) })
  }, error = function(ex) {
    print(ex)
  })
  message(sprintf("*** %s() ... DONE", name))
}

## FIXME:
graphics.off()

message("*** toNnn() ... DONE")
