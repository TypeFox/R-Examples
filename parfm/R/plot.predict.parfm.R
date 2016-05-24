################################################################################
#  Plots the prediction of frailties                                           #
################################################################################
#                                                                              #
#  Computes the prediction of the fraities as                                  #
#                                                                              #
#  Its parameters are                                                          #
#   - x    : the prediction, object of class 'predict.parfm'                   #
#   - sort : how the values must be sorted, either                             #
#            'i' for increasing                                                #
#            'd' for decreasing                                                #
#            or any other value for keeping the current order                  #
#                                                                              #
#   Date: February 02, 2012                                                    #
#   Last modification on: February 23, 2012                                    #
################################################################################

plot.predict.parfm <- function(x, sort="i", 
                               main=NULL,
                               sub=NULL,
                               cex.axis=1,
                               hline=1,
                               ylim=NULL,
                               ...){
  # library('graphics')
  xlab = attr(x, "clustname")
  if (is.null(main)) {
    frailty <- attr(x, "frailty")
    dist <- attr(x, "dist")
    main <- paste(toupper(substr(frailty, 1, 1)), substr(frailty, 2, 100), 
                  " frailty model\nwith ", 
                  toupper(substr(dist, 1, 1)), substr(dist, 2, 100),
                  " baseline", sep="", collapse="")
  }
  if (sort == "i")
    x <- sort(x, decreasing=FALSE)
  else if (sort == "d")
    x <- sort(x, decreasing=TRUE)
  
  if (is.null(ylim))
    ylim <- c(0, max(as.numeric(x)))
  
  plot(1, ty="n", 
       ylab="Predicted frailty value", 
       xlab=xlab,
       main=main,
       sub=sub,
       xaxt="n",
       xlim=c(1, length(as.numeric(x))),
       ylim=ylim)
  abline(h=hline, col="gray")
  points(1:length(as.numeric(x)), as.numeric(x))
  axis(side=1, at=1:length(x), labels=names(x), las=3, cex.axis=cex.axis)
}
