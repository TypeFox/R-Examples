library("R.devices")

message("*** capturePlot() ...")

cat("Default graphics device:\n")
str(getOption("device"))

if (getRversion() >= "3.3.0") {
  g <- capturePlot({
    plot(1:10)
  })

  ## Display
  print(g)

  ## Display with a 2/3 aspect ratio
  toDefault(aspectRatio=2/3, print(g))

  ## Redraw to many output formats
  devEval(c("png", "eps", "pdf"), aspectRatio=2/3, print(g))

} ## if (getRversion() >= "3.3.0")

message("*** capturePlot() ... DONE")
