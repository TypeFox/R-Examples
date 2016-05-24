message("*** dataURI ...")

R.devices::devOptions("*", field="dataURI")
uri <- R.devices::toPNG("foo", tags=c("a", "b"), aspectRatio=0.7, {
  plot(1:10)
})
str(uri)


R.devices::devOptions("*", field="pathname")
png <- R.devices::toPNG("foo", tags=c("a", "b"), aspectRatio=0.7, {
  plot(1:10)
}, force=FALSE)
uri2 <- R.devices::asDataURI(png)
str(uri2)
stopifnot(identical(uri2, uri))

message("*** dataURI ... DONE")
