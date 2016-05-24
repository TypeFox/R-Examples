myFunction = function(nPoints, progressbar=NULL, progresslabel=NULL){
    # Do something that takes a long time
  N = 10
  for(i in 1:N){
    Sys.sleep(0.5)
    if(!missing(progressbar)) progressbar$setFraction(i/N)
    if(!missing(progresslabel)) progresslabel$setText(paste(signif((i)/N*100, 2), "% done"))
  }
  print("Finished")
}

myFunction.dialog = list(show.progress = TRUE, nPoints.integerItem = 10, label = "Number of points to plot")

run.dialog(myFunction)
