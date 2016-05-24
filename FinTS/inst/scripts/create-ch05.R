#### Create the data objects used in chapter 5 
####
####
library(FinTS)
data(TsayFiles)  

##
## 0.  TsayFiles directory 
##

# Adjust the following to
# setwd to 'TsayFiles' 
getwd()
#setwd("..")
#setwd("FinTS")
#setwd("pkg")
#setwd("inst")
#setwd("scripts")
#setwd("TsayFiles")

TsayDir <- "../FinTS/pkg/inst/scripts/TsayFiles/"
dir(TsayDir)

# 0.1.  Tsayfiles$ch05$text 

TsayFiles$ch05$text[1:2,]

ch05 <- with(TsayFiles$ch05, text[text[, 4]=="TRUE", ])
sort(table(ch05[, "data"]))
sort(table(ch05[, "file"]))
# two non-unique pairs:   
# day15.dat and day15.txt
# ibm1to5-dur.dat and ibm1to5-dur.txt

as.vector(ch05[, "data"])

# [1] "ibm"         "ibm9912-tp"  "ibmdurad"    "ibm1to5-dur" "ibm91-ads"  
# [6] "ibm91-adsx"  "day15-ori"   "day15"       "day15"       "ibm1to5-dur"
#[11] "eacd"        "wacd"        "gacd"        "tar-wacd"   

ch05. <- paste(TsayDir, ch05[, "file"], sep="")
(ch05.datNames <- make.names(ch05[, 1]))

##
## 1.  ibm
##     IBM transactions data (11/1/90-1/31/91):
##     The columns are date/time, volume, bid quote, 
##     ask quote, and transaction price
##

readLines(ch05.[1], 4)
ibm. <- readLines(ch05.[1])

# Use chron for days and fractions
# class 'Date' should be whole dayse
library(chron)
ibm0 <- substring(ibm., 1, 10)
ibm0[1:4]
ibm1 <- as.Date(ibm0, "%y%m%d")
ibm1[1:4]
# seconds:  
ibm2 <- substring(ibm., 11, 15)
ibm2[1:4]

# combined:  
ibm3 <- chron(as.numeric(ibm1), as.numeric(ibm2)/(24*3600))
ibm3[1:4]

ibm4 <- read.table(ch05.[1],
         col.names=c("time", "volume", "bid", "ask", "price"))
ibm4[1:4,]
length(ibm3)
length(unique(ibm3))

ibm <- cbind(date.time=ibm3, ibm4[-1])
ibm[1:4,]
##
## 2.  ibm9912-tp
##     IBM transactions data of December 1999.
##     (day, time, price)
##
readLines(ch05.[2], 4)
ibm9912.tp. <- read.table(ch05.[2],
                          col.names=c("day", "seconds", "price"))
ibm9912.tp.[1:2.,]
sapply(ibm9912.tp., range)
#     day seconds   price
#[1,]   1   34205 102.250
#[2,]  31   66507 122.125
ibm9912.tp1 <- with(ibm9912.tp., 
                    chron("11/30/99")+day+(seconds/(24*3600)))
ibm9912.tp1[1:4]
ibm9912.tp <- data.frame(date.time=ibm9912.tp1,
                         price=ibm9912.tp.$price)
ibm9912.tp[1:4, ]

##
## 3.  ibmdurad
##     Adjusted time durations between trades
##     (11/01/90-1/31/91). Positive durations only
##
readLines(ch05.[3], 4)
ibmdurad. <- read.table(ch05.[3],
                       col.names=c("day", "seconds", "adjusted.duration"))
ibmdurad1 <- with(ibmdurad.,
                  chron("11/30/99")+day+(seconds/(24*3600)))
ibmdurad1[1:4]
ibmdurad <- data.frame(date.time=ibmdurad1,
                       adjusted.duration = ibmdurad.$adjusted.duration)
ibmdurad[1:4, ]
##
## 4.  ibm1to5-dur.dat
##     Adjusted durations in (3) for the first 5 trading days
##
ch05.[4]
readLines(ch05.[4], 4)

# Take this from ibmdurad, rather than reading it from the file.

selDur <- with(ibmdurad, ((chron("11/30/99")+1 < date.time)
               & (date.time < chron("11/30/99")+6)
               & (adjusted.duration > 0)) )
ibm1to5.dur <- ibmdurad[selDur, ]
str(ibm1to5.dur)
ibm1to5.dur[1:4,]

# check
ibm1to5.durDat <- scan(ch05.[4])
all.equal(ibm1to5.durDat, ibm1to5.dur$adjusted.duration)
#[1] "Mean relative  difference: 7.657944e-06"
# Comparing the first record from each files
# shows that 'ibmdurad.dat' has 7 significant digits,
# while 'ibm1to5-dur.dat' has only 5.

mean(ibmdurad$adjusted.duration==0)
# 0.1091447

##
## 5.  ibm91-ads.dat
##     Data for Example 5.2:  the ADS file 
##
readLines(ch05.[5], 4)
ibm91.ads <- read.table(ch05.[5],
     col.names=c("A.priceChange", "DirectionOfChg", "SizeInTicks"))
str(ibm91.ads)
first1 <- which(ibm91.ads$SizeInTicks==1)

with(ibm91.ads, sum((A.priceChange!=0)&(SizeInTicks==0)))
# 10

lapply(ibm91.ads, table)
#$V1
#    0     1 
#40057 19718 
#$V2
#   -1     0     1 
# 9835 40057  9883 
#$V3
#    0     1     2     3     4     5     6     7     8     9    10    11    12 
#40067 17384  1555   313   162   106    65    63    25     7     6     1     4 
#   13    14    15    16    21    22    25    29 
#    1     3     1     2     2     2     2     4 

sum(with(ibm91.ads, (A.priceChange!=0)&(SizeInTicks==0)))
# 10

##
## 6.  ibm91-adsx.dat
##     Data for Example 5.2:  The explanatory variables as defined
##
readLines(ch05.[6], 4)
ibm91.adsx <- read.table(ch05.[6], col.names =
      c("Volume.thousands", "time.betw.trades", "bid.ask.spread",
        "A.priceChange", "DirectionOfChg", "SizeInTicks"))

ads.ne.x <- outer(which(ibm91.ads[-59775,1]!= ibm91.adsx[-1, 4]),
                  0:1, "+")
ads.ne <- as.vector(t(ads.ne.x))
  
ibm91.ads[-59775,][ads.ne,]
ibm91.adsx[-1, ][ads.ne,]
table(ibm91.ads[-59775,1], ibm91.adsx[-1, 4])
table(ibm91.ads[-59775,2], ibm91.adsx[-1, 5])
table(ibm91.ads[-59775,3], ibm91.adsx[-1, 6])
table(ibm91.adsx[,3])

##
## 7.  day15-ori.dat
##     Transactions data of IBM stock on November 21, 1990
##     original data
##
readLines(ch05.[7], 4)
day15.ori. <- read.table(ch05.[7])
str(day15.ori.)
# 728 obs on 2 vars
transTm <- (strptime("11/21/1990", "%m/%d/%Y", "ESTSEDT")
            +day15.ori.[[1]])
transTm[1:4]

sum(diff(transTm)==0)
# 60 
day15.ori <- data.frame(transactionTime=transTm,
                        stockPrice=day15.ori.[[2]])
day15.ori[1:2,]
str(day15.ori)
# 728
range(day15.ori$transactionTime)
#[1] "1990-11-21 09:31:17 EST"
#[2] "1990-11-21 16:00:49 EST"

##
## 8.  day15.dat
##     Transactions data of IBM stock on November 21, 1990
##     data for PCD models
##
readLines(ch05.[8], 4)
day15. <- read.table(ch05.[8])
sapply(day15., range)
day15a <- read.table(ch05.[8],
    col.names=c("day", "seconds", "timeBetwPriceChg", "DirectionOfChg", 
      "priceChgTicks", "nTradesWoChg", "multTrans", "dailyCumChg"))
day15a[1:4,]

daySec <- (strptime("11/21/1990", "%m/%d/%Y", "ESTSEDT")
           +day15a$seconds)
daySec[1:11]
sum(diff(daySec)==0)
# 0
day15 <- zoo(day15a[-(1:2)], daySec)
day15[1:11,]
str(day15)
# 195
day15[194:195,]

plot(index(day15), day15[, "dailyCumChg"])
sum(day15[, "nTradesWoChg"])
# 475
# + 194 = 669 NOT 728 ... ?  
sum(day15[, "multTrans"])
# 19 
##
## 9.  day15.txt
##     Transactions data of IBM stock on November 21, 1990
##     data descriptions
##
readLines(ch05.[9])


##
## 10.  ibm1to5-dur.txt
##      The data file used with 
##      RATS programs for estimating duration models:
##   
readLines(ch05.[10], 4)
ibm1to5.durTxt <- scan(ch05.[10])
all.equal(ibm1to5.durDat, ibm1to5.durTxt)
# TRUE

##
## 11.  eacd.rats
##      RATS programs for estimating duration models:
##      EACD model
##

##
## 12.  wacd.rats
##      RATS programs for estimating duration models:
##      WACD model
##

##
## 13.  gacd.rats
##      RATS programs for estimating duration models:
##      GACD model
##

##
## 14.  tar-wacd.rats
##      RATS programs for estimating duration models:
##      Threshold-WACD model
##








##
## 15.  Write data files
##
ch05.rda <- paste(ch05.datNames, "rda", sep=".")

#sel4 <- c(1, 6:8)
sel5 <- 1:8

for(i in sel5)
  save(list=ch05.datNames[i], file=ch05.rda[i])
