plotJGET<-function(J, SHOWONLY=FALSE)
  {
    if(missing(SHOWONLY)) SHOWONLY=FALSE
    ##  Program for plotting the output of JGETseis
##    J  = JGET.seis(fn,kind=2,BIGLONG=FALSE,HEADONLY=FALSE,Iendian=3,PLOT=0)

    GH=prepSEIS(J)
    swig(GH, SHOWONLY=SHOWONLY)

    invisible(GH)
  }



