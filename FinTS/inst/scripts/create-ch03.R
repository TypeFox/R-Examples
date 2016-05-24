#### Create the data objects used in chapter 3 
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

# 0.1.  Tsayfiles$ch03$text 

TsayFiles$ch03$text[1:2,]

ch03 <- with(TsayFiles$ch03, text[text[, 4]=="TRUE", ])
sort(table(ch03[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch03[, "data"])

#[1] "m-intc7303"      "exch-perc"       "sp500"           "m-ibm2697"      
#[5] "d-ibmvwewsp6203" "m-ibmspln"       "m-ibmsplnsu"     "d-sp8099"       

ch03. <- paste(TsayDir, ch03[, "file"], sep="")
(ch03.datNames <- make.names(ch03[, 1]))

##
## 1.  m-intc7303
##     Monthly simple returns of Intel stock
##

readLines(ch03.[1], 4)
tst.intc7303 <- read.zoo(ch03.[1], format="%Y%m%d",
      col.names=c("Date", "Intel"))
data(m.intc7303)
all.equal(tst.intc7303, m.intc7303)
# TRUE 
##
## 2.  exch-perc
##     10-minute FX log returns (Mark-Dollar)
##
readLines(ch03.[2], 4)
exch.perc <- scan(ch03.[2])
# read 2497 items;  book says 2488 :-(
sum(is.na(exch.perc))
# 0 
plot(exch.perc, type="l")

##
## 3.  sp500
##     Monthly excess returns of the S&P 500 index
##
readLines(ch03.[3], 4)
sp500a <- scan(ch03.[3])
# read 792 items starting from 1926 
# (91-47)*4 = 864
sp500 <- zooreg(sp500a, 1926, freq=12)
sp500[1:4]

##
## 4.  m-ibm2697
##     Monthly returns of IBM stock
##
readLines(ch03.[4], 4)
tst.ibm2697a <- scan(ch03.[4])
# read 864 items
# 97-25 = 72 * 12 = 864 

tst.ibm2697 <- zooreg(tst.ibm2697a, start=1926, freq=12)
                
data(m.ibm2697)
all.equal(tst.ibm2697, m.ibm2697)
# TRUE 
##
## 5.  d-ibmvwewsp6203
##     Daily returns of IBM stock, VW, EW, and SP500
##
readLines(ch03.[5], 4)
tst.ibmvwewsp6203 <- read.zoo(ch03.[5], format="%Y%m%d",
      col.names=c("Date", "IBM", "VW", "EW", "SP") )
tst.ibmvwewsp6203[1:4,]
range(index(tst.ibmvwewsp6203))
data(d.ibmvwewsp6203)
range(index(d.ibmvwewsp6203))

all.equal(tst.ibmvwewsp6203, d.ibmvwewsp6203)
# TRUE 
##
## 6.  m-ibmspln
##     Monthly log returns of IBM stock and S&P 500 index
##
readLines(ch03.[6], 4)
m.ibmspln. <- read.table(ch03.[6], col.names=c("IBM", "SP"))
dim(m.ibmspln.)
#  888  2
m.ibmspln.[1:4,]
m.ibmspln <- zooreg(m.ibmspln., start=1926, freq=12)
m.ibmspln[1:4,]
str(m.ibmspln)
##
## 7.  m-ibmsplnsu
##     Data for Example 3.4
##
readLines(ch03.[7], 4)
m.ibmsplnsu. <- read.table(ch03.[7],
                col.names=c("IBM", "SP", "summer"))
                                      
m.ibmsplnsu <- zooreg(m.ibmsplnsu., start=1926, freq=12)
m.ibmsplnsu[1:4,]
str(m.ibmsplnsu)

table(m.ibmsplnsu[, "summer"])
plot(m.ibmsplnsu[, "summer"])
plot(m.ibmsplnsu[1:33, "summer"])

##
## 8.  d.sp8099
##     Daily returns of S&P 500 index:
##
readLines(ch03.[8], 4)
d.sp8099 <- read.zoo(ch03.[8], format="%Y%m%d",
                   col.names=c("Date", "SP") )
d.sp8099[1:4,]
range(index(d.sp8099))

##
## 9.  Write the data files
##

ch03.rda <- paste(ch03.datNames, "rda", sep=".")

sel3 <- c(2:3, 6:8)

for(i in sel3)
  save(list=ch03.datNames[i], file=ch03.rda[i])
