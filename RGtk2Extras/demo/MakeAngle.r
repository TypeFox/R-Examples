  # example 3: an angle view displayer, from Graham Williams
  makeAngle <-function(unit="degrees",
                       degrees=36,
                       radians=0.63,
                       turns=0.1,
                       sine=c("Do not show"),
                       delay=10)
  {
    if (degrees != 36)
    {
      radians <- round(degrees/180*pi, 2)
      turns <- round(degrees/360, 3)
    }
    else if (radians != 0.63)
    {
      degrees <- round(radians*180/pi)
      turns <- round(radians/(2*pi),2)
    }
    else if (turns != 0.1)
    {
      degrees <- round(turns*360)
      radians <- round(turns*2*pi)
    }

    numsteps <- 50
    angle <- c(degrees=degrees, radians=radians, turns=turns)[unit]
    if(abs(angle)>1)meas=paste(unit,"s", sep="") else meas=unit
    opar <- par(mar=rep(2,1,4), mgp=c(2,.5,0), las=1)
    plot(c(-1,1), c(-1,1), type="n", bty="n", xlab="", ylab="", axes=F)
    axis(1, pos=0, col="gray", at=c(-1, 1))
    axis(2, pos=0, col="gray", at=c(-1, 1), srt=-90)
    symbols(x=0, y=0, circles=1, add=T, inches=F)
    lines(c(0,1), c(0,0), col="red", lwd=2)
    theta <- radians
    bits <- seq(from=0, to=theta, length=numsteps)
    for(i in 2:length(bits))
    {
      ci0 <- cos(bits[i-1])
      si0 <- sin(bits[i-1])
      ci1 <- cos(bits[i])
      si1 <- sin(bits[i])
      segments(ci0, si0, ci1, si1, col="red", lwd=2)
      Sys.sleep(0.001*delay)
      if(i==length(bits)) lines(c(0,ci1), c(0, si1), col="red", lwd=2)
    }
    for(i in 2:numsteps){
      ci0 <- cos(bits[i-1])
      si0 <- sin(bits[i-1])
      ci1 <- cos(bits[i])
      si1 <- sin(bits[i])
      segments(0.25*ci0, 0.25*si0, 0.25*ci1, 0.25*si1, col="gray", lwd=2)
    }
    arrows(0.25*cos(bits[numsteps-1]), 0.25*sin(bits[numsteps-1]),
           0.25*cos(bits[numsteps]), 0.25*sin(bits[numsteps]))
    srot <- theta*90/pi
    end <- 0
    if(cos(theta/2)<0){
      srot <- srot+180
      end=1
    }
    text(0.25*cos(theta/2), .25*sin(theta/2), adj=end,
         substitute(theta==z*phantom(1)*lab, list(z=angle, lab=meas)),
         srt=srot)
    if(sine=="Show"){
      ch <- strheight("0")
      ctheta <- cos(theta)
      stheta <- sin(theta)
      srot <- theta*180/pi
      if(ctheta<0) srot <- srot+180
        text(0.5*cos(theta)-0.75*ch*sin(theta), 0.5*sin(theta)+0.75*ch*cos(theta),
             "Hypotenuse = 1", adj=0.5, srt=srot)
      lines(c(cos(theta), cos(theta)), c(0,sin(theta)))
      lines(c(0,cos(theta)), c(0,0))
      if(sin(theta)!=0)
        text(cos(theta)-0.85*strwidth("0"),0.5*sin(theta),
             substitute(sin*phantom(".")*theta == si,
                        list(si=round(stheta,3))),
             adj=0.5, srt=90)
      if(cos(theta)!=0)
        text(0.5*cos(theta),-0.85*strheight("0"),
             substitute(cos*phantom(".")*theta  == ci,
                        list(ci=round(ctheta,3))), adj=0.5)
    }
    return(c(degrees, radians, turns))
  }
  
  update_ranges = function(dummy, degrees, radians, turns, user.data){
    theUnits <- c(degrees=1, radians=180/pi, turns=360)
    theWidgets <- c(degrees=degrees, radians=radians, turns=turns)
    theVals <- sapply(theWidgets, get.value) * theUnits
      # Turn propagate off here to avoid infinite loops
    for(nn in setdiff(names(theWidgets), user.data))
      set.value(theWidgets[[nn]], theVals[user.data]/theUnits[nn], propagate=FALSE)
    return(FALSE)
  }

   makeAngle.dialog <- list(label="Demonstrate Angles",
     keep.open = TRUE,
     unit.radiobuttonItem=c(value="degrees", "radians", "turns"), item.labels = c("Degrees", "Radians", "Turns"), label = "Unit",
     tooltip="Use these units on the graph.",
     signal = c("default", function(unit, degrees, radians, turns){
        widgets <- c(degrees=degrees, radians=radians, turns=turns)
        sapply(widgets, function(x) x$setSensitive(FALSE))
        widgets[[get.value(unit)]]$setSensitive(TRUE)
     }, "degrees", "radians", "turns"),

     degrees.rangeItem=c(value=36, from=-360, to=360, by=3),
       signal = list("default", update_ranges, "degrees", "radians", "turns", user.data="degrees"),
     label="Degrees",
     tooltip="Angle in degrees.",

     radians.rangeItem=c(value=0.63, from=-6.283, to=6.283, by=0.05236),
       signal = list("default", update_ranges, "degrees", "radians", "turns", user.data="radians"),
     label="Radians",
     tooltip="Angle in radians.",

     turns.rangeItem=c(value=0.1, from=-1, to=1, by=1/120),
       signal = list("default", update_ranges, "degrees", "radians", "turns", user.data="turns"),
     label="Turns",
     tooltip="Anti-clockwise turns.",

     sine.choiceItem=c(value="Do not show", "Show"),
     label="Show sine",
     tooltip="Show sine and cosine of angle.",

     delay.integerItem=3,
     tooltip = "Delay time (ms)",
     label="Delay")

run.dialog(makeAngle)

  ## example 4:
  ##