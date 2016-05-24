#### Create the data objects used in chapter 1 
####
####
library(FinTS)
data(TsayFiles)  

#*** Thanks to Gabor Grothendieck and Diethelm Wuertz
#*** for help with how to process dates and times 

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

# 0.1.  Tsayfiles$ch01$text 


TsayFiles$ch01$text[1:2,]
ch01 <- with(TsayFiles$ch01, text[text[, 4]=="TRUE", ])
sort(table(ch01[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch01[, "data"])

# [1] "d-ibmvwewsp6203" "d-intc7303"      "d-3m6203"        "d-msft8603"     
# [5] "d-c8603"         "m-ibmvwewsp2603" "m-intc7303"      "m-3m4603"       
# [9] "m-msft8603"      "m-c8603"         "m-gs10"          "m-gs1"          
#[13] "d-fxjp00"        "m-fama-bond5203" "m-gs3"           "m-gs5"          
#[17] "w-tb3ms"         "w-tb6ms"

ch01. <- paste(TsayDir, ch01[, "file"], sep="")

##
## 1.  d-ibmvwewsp6203
##     Daily simple returns of IBM, VW, EW, SP (7/3/62-12/31/03):
##      (Format: date, IBM, VW, EW & SP)
##

d.ibmvwewsp6203 <- read.zoo(ch01.[1], format="%Y%m%d",
           col.names=c("date", "IBM", "VW", "EW", "SP") )
d.ibmvwewsp6203[1:2,]
class(index(d.ibmvwewsp6203))
# Date 
##
## 2.  d-intc7303
##     Daily simple returns of Intel stock (12/15/72-12/31/03) 
##
d.intc7303 <- read.zoo(ch01.[2], format="%Y%m%d",
           col.names=c("date", "Intel") )
d.intc7303[1:2, ]
index(d.intc7303)[1:2]
class(index(d.intc7303))
# Date
##
## 3.  d-3m6203
##     Daily simple returns of 3M stock
##
d.3m6203 <- read.zoo(ch01.[3], format="%Y%m%d",
           col.names=c("date", "3M") )
d.3m6203[1:2, ]

##
## 4.  d-msft8603
##     Daily simple returns of Microsoft stock
##
d.msft8603 <- read.zoo(ch01.[4], format="%Y%m%d",
           col.names=c("date", "msft") )
d.msft8603[1:2, ]

##
## 5.  d-c8603
##     Daily simple returns of Citi-group stock
##
d.c8603 <- read.zoo(ch01.[5], format="%Y%m%d",
           col.names=c("date", "citi") )  
d.c8603[1:2]

##
## 6.  m-ibmvwewsp2603
##     Monthly simple returns of IBM, VW, EW, SP (1/26-12/03)
##     (Format: date, IBM, VW, EW, & SP) 
m.ibmvwewsp2603 <- read.yearmon(ch01.[6], format="%Y%m%d",
           col.names=c("date", "IBM", "VW", "EW", "SP") )  
m.ibmvwewsp2603[1:2,]
index(m.ibmvwewsp2603)[1:2]

##
## 7.  m-intc7303
##     Monthly simple returns of Intel stock
##
m.intc7303 <- read.yearmon(ch01.[7], format="%Y%m%d",
           col.names=c("date", "Intel") )
m.intc7303[1:2, ]
index(m.intc7303)[1:2]
##
## 8.  m-3m4603
##     Monthly simple returns of 3M stock
##
m.3m4603 <- read.yearmon(ch01.[8], format="%Y%m%d",
           col.names=c("date", "3M") ) 
m.3m4603[1:2,]
##
## 9.  m-msft8603
##     Monthly simple returns of Microsoft stock
##
m.msft8603 <- read.yearmon(ch01.[9], format="%Y%m%d",
           col.names=c("date", "msft") ) 
m.msft8603[1:2,]
##
## 10.  m-c8603
##      Monthly simple returns of Citi-group stock
##
m.c8603 <- read.yearmon(ch01.[10], format="%Y%m%d",
           col.names=c("date", "citi") )
m.c8603[1:2, ]
##
## 11.  m-gs10 & m-gs1
##      Monthly 10-yr  and 1-yr Treasury constant maturity rates
##      (4/53-3/04):  (Format: year, month, date, rate): 
##
#  Nonstandard date format makes reading difficult 
ch01[11:12,]
readLines(ch01.[11], 4)
#m.gs10 <- read.zoo(ch01.[11], "%Y %m %d",col.names=c("date","gs10"))
# Error:  More columns than column names

m.gs10a <- readLines(ch01.[11])
m.gs10a[1:2]
m.gs10t <- as.Date(m.gs10a, "%Y %m %d") 
m.gs10t[1:4]

m.gs10 <- zoo(as.numeric(substring(m.gs10a, 11)),
              as.yearmon2(m.gs10t) )
                
m.gs10[1:2,]
index(m.gs10)[1:2]

m.gs1a <- readLines(ch01.[12])
m.gs1a[1:2]
m.gs1 <- zoo(as.numeric(substring(m.gs1a, 11)),
             as.yearmon2(as.Date(m.gs1a, "%Y %m %d")))
m.gs1[1:2,]
index(m.gs1)[1:2]

##
## 12.  d-fxjp00
##      Daily exchange rate between U.S. dollar and Japanese yen
##      (Format: ddmmyy, fx)
##
ch01[13,]
readLines(ch01.[13], 4)
#d.fxjp00 <- read.zoo(ch01.[13], format="%d%m%y",
#           col.names=c("date", "fxjp") )
# ERROR:  index contains NAs
# Parse manually
d.fxjp00a <- readLines(ch01.[13])
d.fxjp00a[1:4]
# Insert sep characters
# because I can't seem to get strptime to work properly
# otherwise with this 2-digit year format.
d.fxjp00b <- paste(substring(d.fxjp00a, 1, 2),
                   substring(d.fxjp00a, 3), sep="-")
d.fxjp00t <- as.Date(d.fxjp00b, "%d-%m%y")
d.fxjp00t[1:4]
class(d.fxjp00t)

d.fxjp00 <- zoo(as.numeric(substring(d.fxjp00a, 7)),
                 d.fxjp00t ) 
d.fxjp00[1:4] # good 

##
## 13.  m-fama-bond5203
##      Monthly bond returns (1-12m, 24-36m, 48-60m, 61-120m)
##      (Format: date, bond returns) 
##
ch01[14, ]
readLines(ch01.[14], 4)
m.fama.bond5203 <- read.yearmon(ch01.[14], format="%Y%m%d",
      col.names=c("date", "m1.12", "m24.36", "m48.60", "m61.120") )
m.fama.bond5203[1:2, ]
##
## 14.  m-gs3 & m-gs5 
##      Monthly 3-yr and 5-yr Treasury constant maturity rates
##
ch01[15:16,]
readLines(ch01.[15], 4)
# Problem with another non-standard date format
m.gs3a <- readLines(ch01.[15])
m.gs3d <- as.Date(m.gs3a, "%Y %m %d")
m.gs3d[1:4]

m.gs3 <- zoo(as.numeric(substring(m.gs3a, 11)),
             as.yearmon2(m.gs3d))
m.gs3[1:4]

m.gs5a <- readLines(ch01.[16])
m.gs5a[1:4]
m.gs5d <- as.Date(m.gs5a, "%Y %m %d")
m.gs5d[1:4]
m.gs5 <- zoo(as.numeric(substring(m.gs5a, 11)),
             as.yearmon2(m.gs5d))
m.gs5[1:4]

##
## 15.  w-tb3ms & w-tb6ms 
##      Weekly Treasury Bill rates - 3 & 6 months
##
ch01[17:18, ]
w.tb3ms <- read.zoo(ch01.[17], format="%Y%m%d",
           col.names=c("date", "tb3ms") )
w.tb3ms[1:2, ]
w.tb6ms <- read.zoo(ch01.[18], format="%Y%m%d",
           col.names=c("date", "tb3ms") ) 
w.tb6ms[1:2,]

##
## 16.  Write the data file
##

ch01.datNames <- make.names(ch01[, 1])
ch01.rda <- paste(ch01.datNames, "rda", sep=".")
#save(list=ch01.datNames, file="ch01.rda")

nObj1 <- length(ch01.datNames)
for(i in 1:nObj1)
  save(list=ch01.datNames[i],
       file=ch01.rda[i])

#####################
# used to fix m.* after read.yearmon and as.yearmon2 were written:  
#m.obj <- c(6:12, 14:16)
#for(i in m.obj)save(list=ch01.datNames[i], file=ch01.rda[i])
