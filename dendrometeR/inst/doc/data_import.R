## ---- include = FALSE----------------------------------------------------
library(dendrometeR)

## ---- echo = FALSE, results = 'asis'-------------------------------------
data(dmCDraw)
knitr::kable(head(dmCDraw, 5))

## ---- eval = FALSE-------------------------------------------------------
#  dm.data <- data.frame(timestamp = ISOdate(year = dmCDraw$year, month = dmCDraw$Month,
#  day = dmCDraw$Day, hour = dmCDraw$H, min = dmCDraw$M))

## ---- eval = FALSE-------------------------------------------------------
#  dm.data$dendro <- dmCDraw$dendro

## ---- eval = FALSE-------------------------------------------------------
#  rownames(dm.data) <- dm.data$timestamp
#  dm.data$timestamp <- NULL

## ---- echo = FALSE, results = 'asis'-------------------------------------
data(dmCD)
knitr::kable(head(dmCD, 5))

## ---- echo = FALSE, results = 'asis'-------------------------------------
data(dmHSraw)
knitr::kable(head(dmHSraw, 5))

## ---- eval = FALSE-------------------------------------------------------
#  dm.data <- data.frame(timestamp = as.POSIXct(strptime(dmHSraw$TIMESTAMP, '%Y-%m-%d %H:%M:%S'), tz = "GMT"))

## ---- eval = FALSE-------------------------------------------------------
#  dm.data$dBUP2 <- dmHSraw$dBUP2

## ---- eval = FALSE-------------------------------------------------------
#  rownames(dm.data) <- dm.data$timestamp
#  dm.data$timestamp <- NULL

## ---- echo = FALSE, results = 'asis'-------------------------------------
data(dmEDraw)
knitr::kable(head(dmEDraw, 5))

## ---- eval = FALSE-------------------------------------------------------
#  dm.data <- data.frame(timestamp = as.POSIXct(strptime(dmEDraw$TIMESTAMP, '%Y-%m-%d %H:%M:%S'), tz = "GMT"))

## ---- eval = FALSE-------------------------------------------------------
#  # option 1: select series by typing column names:
#  dm.data[,c("Beech03","Beech04")] <- dmEDraw[,c(2,3)]
#  
#  # option 2: select series from the character vector produced by names:
#  dm.data[,names(dmEDraw)[c(2,3)]] <- dmEDraw[,c(2,3)]

## ---- eval = FALSE-------------------------------------------------------
#  dm.data[,names(dmEDraw)[2:3]] <- dmEDraw[,2:3]

## ---- eval = FALSE-------------------------------------------------------
#  rownames(dm.data) <- dm.data$timestamp
#  dm.data$timestamp <- NULL

## ------------------------------------------------------------------------
is.dendro(dmCDraw)
is.dendro(dmCD)

## ------------------------------------------------------------------------
dendro.resolution(dmCD)
dendro.resolution(dmCD, unts = "hours")

## ------------------------------------------------------------------------
TRUE %in% is.na(dmCD)
data(dmED)
TRUE %in% is.na(dmED)

