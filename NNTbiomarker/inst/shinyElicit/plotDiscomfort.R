
plotDiscomfort = function(
    NNTlower = 7,
    NNTupper = 16,
    NNTpos = 1,
    NNTneg=NNTupper*1.5,
    drawAxes=T,
    drawNNT=T,
    drawPosNeg=T){

  labelCex = 2.0
  triangleUpPch = 24
  triangleDownPch = 25
  numbersY = 0.5
  YtriangleDown = 1.5
  YnntLowerUpper = 2.5
  barY = 3; barHeight = 1
  YtriangleUp = 4
  YnntPosNeg = 5
  rightSideLimit = max(NNTneg * 1.15, NNTneg + 3, NNTupper * 1.15, NNTupper + 3)
  par( mar=c(0.01,0.01,0.01,0.01) #mgp=c(0,0,0)
       #oma=c(0,0,0,0) # does nothing
       #omd=c(0,0,0,0) # fig margins too wide
       )
  YlimUpper = ifelse(drawPosNeg, 7, 4)
  plot(c(1,rightSideLimit), c(0, YlimUpper), pch="",
       axes=drawAxes, xlab="", ylab=""
       )
  actBarLeft = 1
  actBarRight = NNTlower
  actBarX = actBarLeft + (actBarRight - actBarLeft)/2
  actBarWidth = actBarRight - actBarLeft
  symbols(actBarX, barY, inches=F,
          rectangles=matrix(c(actBarWidth, barHeight), nrow=1), bg="green", add=T)
  graphics::text(actBarX, barY, "Act!", cex=labelCex)

  discomfortBarLeft = NNTlower
  discomfortBarRight = NNTupper
  discomfortBarX = discomfortBarLeft + (discomfortBarRight - discomfortBarLeft)/2
  discomfortBarWidth = discomfortBarRight - discomfortBarLeft
  symbols(discomfortBarX, barY, inches=F,
          rectangles=matrix(c(discomfortBarWidth, barHeight), nrow=1),
          bg="red", add=T)
  graphics::text(discomfortBarX, barY, "discomfort range", cex=labelCex)

  waitBarLeft = NNTupper
  waitBarRight = rightSideLimit
  waitBarX = waitBarLeft + (waitBarRight - waitBarLeft)/2
  waitBarWidth = waitBarRight - waitBarLeft
  symbols(waitBarX, barY, inches=F,
          rectangles=matrix(c(waitBarWidth, barHeight), nrow=1),
          bg="green", add=T)
  graphics::text(waitBarX, barY, "Wait!", cex=  labelCex)

  numberGap = ceiling(rightSideLimit/25)
  numberSequence = seq(1, rightSideLimit, numberGap)
  graphics::text(x = numberSequence,y = numbersY,
       labels = numberSequence, cex=labelCex,
       col=c("black", "red")
       [1+numberSequence %between% c(NNTlower, NNTupper)]
  )
#   symbols(c(NNTupper, NNTlower), c(numbersY,numbersY),
#           circles=c(0.4,0.4),
#           inches=F, add=T, fg="red", xpd=F)
  points(NNTlower, YtriangleDown, cex=3, pch=triangleDownPch, bg="red")
  points(NNTupper, YtriangleDown, cex=3, pch=triangleDownPch, bg="red")
  graphics::text(discomfortBarLeft, YnntLowerUpper,
       "NNTlower", pos=1, cex=2, col="red")
  graphics::text(discomfortBarRight, YnntLowerUpper,
       "NNTupper", pos=1, cex=2, col="red")
  if(drawNNT) {
    NNT = 1/input$prevalence
    YprevCircle = YnntLowerUpper
    prevCirclePch = 19  # solid circle; filled circle=21
    graphics::text(NNT, YnntLowerUpper, pos=1, "NNT", cex=2)
    points(NNT, YtriangleDown, cex=3, pch=triangleDownPch, bg="black")
    #points(NNT, YprevCircle, cex=5, pch=prevCirclePch, bg="black")
    #points(NNT, YprevCircle, cex=3, pch="P", col="white")
  }
  if(drawPosNeg) {
    points(NNTpos, YtriangleUp, cex=3, pch=triangleUpPch, bg="darkgreen")
    points(NNTneg, YtriangleUp, cex=3, pch=triangleUpPch, bg="darkgreen")
    graphics::text(NNTpos, YnntPosNeg,
         "NNTpos", pos=3, cex=2, col="darkgreen")
    graphics::text(NNTneg, YnntPosNeg,
         "NNTneg", pos=3, cex=2, col="darkgreen")
  }
}

#plotDiscomfort()
