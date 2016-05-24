## start up shiny here
if(!require(shiny)) {
  message("The `shiny` package is required to view these demos")
  stop()
}

demos <- read.csv(textConnection("
Description, fname
Influence of a point on the mean,mean-influence
Histogram and bin selection,hist-bins
Kernel and bandwith selection for a density plot,bw-selection
"))


out <- menu(demos$Description,
            graphics=FALSE,
            title="Demos aninmated through the `shiny` package.

Please select one of the following.
These require that the `shiny` package be installed prior to usage.")

if(out > 0) {
  app <- demos$fname[out]
  f <- system.file("shiny", app, package="UsingR")
  runApp(f)
}
