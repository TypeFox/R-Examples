library(shiny)
library(shotGroups)

unitsDst    <- c("m"="1",  "yard"="2", "feet"="3")
unitsDstInv <- c("1"="m",  "2"="yd", "3"="ft")

unitsXY     <- c("cm"="1", "mm"="2", "inch"="3")
unitsXYInv  <- c("1"="cm", "2"="mm", "3"="in")

unitsAbs    <- c("m"="1", "cm"="2", "mm"="3", "yard"="4", "feet"="5", "inch"="6")
unitsAbsInv <- c("1"="m", "2"="cm", "3"="mm", "4"="yd", "5"="ft", "6"="in")

unitsAng    <- c("degree"="1", "MOA"="2", "SMOA"="3", "radian"="4", "milliradian"="5", "NATO mil"="6")
unitsAngInv <- c("1"="deg",    "2"="MOA", "3"="SMOA", "4"="rad",    "5"="mrad",        "6"="mil")
