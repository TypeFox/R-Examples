#################################################
### code chunk number 136: BasicComponents.Rnw:769-779
###################################################
window <- gtkWindow(); window$setTitle("Progress bar example")
progress_bar <- gtkProgressBar()
window$add(progress_bar)
#
progress_bar$setText("Please be patient...")
for(i in 1:100) {
  progress_bar$setFraction(i/100)
  Sys.sleep(0.05) ## replace with a step in the process
}
progress_bar$setText("All done.")


###################################################
### code chunk number 137: gtk-widget-progress-pulse
###################################################
progress_bar$pulse()

