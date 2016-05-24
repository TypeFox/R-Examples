### tests.write.fwf.R
###------------------------------------------------------------------------
### What: Tests for write.fwf
### $Id: tests.write.fwf.R 1300 2008-08-05 11:47:18Z ggorjan $
### Time-stamp: <2008-08-05 12:25:05 ggorjan>
###------------------------------------------------------------------------

library(gdata)

## --- Test data ---

num <- round(c(733070.345678, 1214213.78765456, 553823.798765678,
               1085022.8876545678,  571063.88765456, 606718.3876545678,
               1053686.6, 971024.187656, 631193.398765456, 879431.1),
             digits=3)

testData <- data.frame(num1=c(1:10, NA),
                       num2=c(NA, seq(from=1, to=5.5, by=0.5)),
                       num3=c(NA, num),
                       int1=c(as.integer(1:4), NA, as.integer(4:9)),
                       fac1=factor(c(NA, letters[1:9], "hjh")),
                       fac2=factor(c(letters[6:15], NA)),
                       cha1=c(letters[17:26], NA),
                       cha2=c(NA, "longer", letters[25:17]),
                       stringsAsFactors=FALSE)
levels(testData$fac1) <- c(levels(testData$fac1), "unusedLevel")
testData$Date <- as.Date("1900-1-1")
testData$Date[2] <- NA
testData$POSIXt <- as.POSIXct(strptime("1900-1-1 01:01:01", format="%Y-%m-%d %H:%M:%S"))
testData$POSIXt[5] <- NA

## --- Tests ---

## Default
write.fwf(testData)

## NA should be -
write.fwf(x=testData, na="-")
## NA should be -NA-
write.fwf(x=testData, na="-NA-")

## Some other separator than space
write.fwf(testData[, 1:4], sep="-mySep-")

## Justify right for character columns
write.fwf(testData, justify="right")

## Justify right for character columns, track na values
write.fwf(testData, justify="right", na="!")

## With quotes
write.fwf(testData, quote=TRUE)

## With rownames
write.fwf(testData, rownames=TRUE)

## Without colnames
write.fwf(testData, colnames=FALSE)

## Without colnames, with rownames
write.fwf(testData, colnames=FALSE, rownames=TRUE)

## With rownames and colnames and rowCol
write.fwf(testData, rownames=TRUE, rowCol="HI!")

## formatInfo is tested with unit tests

###------------------------------------------------------------------------
### tests.write.fwf.R ends
