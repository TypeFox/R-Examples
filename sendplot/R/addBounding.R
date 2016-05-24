#
# add bounding points to a given figure 
#

addBounding <- function(Splot,
                        figure,
                        bb.clr = "blue",
                        bb.cex = 2,
                        boundFileName = "SplotDot",
                        dir="./"
                        ){



  #
  # now make image with bouding points for given figure
  #

  # orig.plt.extras = Splot$plt.extras

  # create call that will be temporarily used in Splot plt.extras
  #   by accessing specified figures bounding limits
  addPts = paste("points(", Splot$plot.lims$xmins[figure], ",", Splot$plot.lims$ymaxs[figure], ",pch=22,bg='", bb.clr, "', col='", bb.clr, "',cex=", bb.cex, ");points(", Splot$plot.lims$xmaxs[figure], ",", Splot$plot.lims$ymins[figure], ",pch=22,bg='", bb.clr, "', col='", bb.clr, "',cex=", bb.cex, ")", sep="")
  
  # cat("addBounding: addPts=", addPts, fill=T)
      
  # format plt.extras properly and add additional call 
  pel = length(Splot$plot.extras[[figure]])
  # if figure plt.extras is a single entry
  if(pel == 1){
    # if it is NA can add to figure directly
    if(is.na(Splot$plot.extras[[figure]][1])){
      Splot$plot.extras[[figure]][1] = addPts
    }else{
      # if it already a single entry, make sure the object is a list and add new plt.call
      if( class(Splot$plot.extras[[figure]][1]) != "list")  class(Splot$plot.extras[[figure]][1]) = "list"
      Splot$plot.extras[[figure]][2] = addPts
    }
  # if figure already has multiple calls        
  }else{
    # check that figure object is a list and add new plt.call
    if(class(Splot$plot.extras[[figure]]) != "list") class(Splot$plot.extras[[figure]]) = "list"      
    Splot$plot.extras[[figure]][(pel+1)] = addPts
  }

  # make output with figure having bound points 
  makeSplot(Splot, fname.root=boundFileName, dir=dir, makeInteractive=FALSE, overwriteSourcePlot=c("png", "tiff"), returnObj=TRUE)
  

}
