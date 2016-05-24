#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          viewData.r
#                ---------------
# Author:        Joern Schulz
# Stand:         15.08.2006
#
#########################################################################
# 
#
viewData <- function(Data, Title=NULL, curWindow=TRUE, colors = gray((0:255)/255), s=1)
{
# ========================================================================
#
# Description:
# This function print 1 up to 6 image in a new window on the display or in
# an current window curWindow. The function used the R function 'image'
# to display the image. Note the condition of this function.
#
# DataList  Has to be a list, that contains alternating the the image and
#           the title of the image. One, two, op to six image with title
#           are possible.
# curWindow If curWindow=NULL then a new window is open on the display.
#           Otherwise the current window will be used.
# colors    The colormap, that will be used.
#
# ========================================================================

  
  if(!(is.matrix(Data)) && !(is.list(Data)))
    stop("'Data' has to be of type 'matrix' or 'list'.")
  if(is.matrix(Data)){
    Data <- list(Data)
    lData <- 1
  } else
    lData <- length(Data)
  
  if(is.null(Title))
    Title <- as.list(paste(as.character(1:lData),".Image",sep=""))
  
  if(!(is.character(Title)) && !(is.list(Title)))
    stop("'Title' has to be of type 'character' or 'list'.")
  if(is.character(Title)){
    Title <- list(Title)
    lTitle <- 1
  } else
    lTitle <- length(Title)
  
  if(lData != lTitle){
    cat("The length of list 'Data' is different to 'Title'. No title is used.\n")
    Title <- as.list(paste(as.character(1:lData),".Image",sep=""))
  }
  
  for(i in 1:(lData)){
    if( !(is.matrix(Data[[i]])) ){
      cat("'Data[[", i, "]]' is not a matrix. \n",sep="")
      stop()
    }
    if( !(is.character(Title[[i]])) ){
      cat("'Title[[", 2*i-1, "]]' is not a character. \n",sep="")
      stop()
    }
  }


 # ========================================================================
 # visualization of the data
 #

 ########################################
 # differents Cases
 #
  if(lData == 1){
    width=s*7
    height=s*7
    mfrow=c(1,1)
  } else if(lData == 2){
    width=s*10
    height=s*5
    mfrow=c(1,2)
  } else if(lData == 3){
    width=s*12
    height=s*4.2
    mfrow=c(1,3)
  } else if(lData == 4){
    width=s*8
    height=s*8
    mfrow=c(2,2)
  } else if(lData == 5){
    width=s*12
    height=s*8
    mfrow=c(2,3)
  } else if(lData == 6){
    width=s*12
    height=s*8
    mfrow=c(2,3)
  } else if(lData == 7){
    width=s*13
    height=s*7
    mfrow=c(2,4)
  } else if(lData == 8){
    width=s*13
    height=s*7
    mfrow=c(2,4)
  } else
    stop("Only 1 up to 8 images are supported with this function.")

 ########################################
 # Setting paramater and print image
 #
  if (curWindow == TRUE){
    curWindow <- dev.cur()
    if(curWindow == 1) curWindow<-FALSE
    else dev.set(curWindow)
  } else if (curWindow != FALSE)
    dev.set(curWindow)
   
  if (curWindow == FALSE) {
    if (.Platform$OS.type == "unix")
      X11(width=width, height=height)
    else if(.Platform$OS.type == "windows")
      windows(width=width, height=height)
    else 
      stop("The programm isn't supported on platform-type: '", .Platform$OS.type,"'.")
  }

 # las = Axenbeschriftung horizontal, oma = außerer Rand, mar = innerer Rand
  par(mfrow=mfrow, las=1, oma = c(2,2,2,2)+0.1, mar=c(1,1,3,0.25), mgp=c(3,1,0))
  par(adj=0.5) # Ausrichtung des Text in der Grafik linkbündig=0, mitte=0.5, r=1
  par(plt=c(0.1,0.9,0.1,0.9)) # bereich setzen wo grafik gestzt werden soll


  for(i in 1:(lData)){
    image(Data[[i]], col=colors, axes=FALSE)
    title(Title[[i]])
  }

  
  invisible(dev.cur())
}
