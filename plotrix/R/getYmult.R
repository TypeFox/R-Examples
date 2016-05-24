getYmult<-function() {
 if(dev.cur() == 1) {
  warning("No graphics device open.")
  ymult<-1
 }
 else {
  # get the plot aspect ratio
  xyasp<-par("pin")
  # get the plot coordinate ratio
  xycr<-diff(par("usr"))[c(1,3)]
  ymult<-xyasp[1]/xyasp[2]*xycr[2]/xycr[1]
 }
 return(ymult)
}
