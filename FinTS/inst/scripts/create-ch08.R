#### Create the data objects used in chapter 8 
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

# 0.1.  Tsayfiles$ch08$text 

TsayFiles$ch08$text[1:2,]

ch08 <- with(TsayFiles$ch08, text[text[, 4]=="TRUE", ])
sort(table(ch08[, "data"]))
# m-ibmspln appears twice
TsayFiles$ch08$text[, c(1:2, 4)]
# m-ibmspln.txt is referenced in the web page
# m-ibmspln.dat is not visible ... 

as.vector(ch08[, "data"])
# [1] "m-ibmsp2699"  "m-ibmspln"    "m-ibmspln"    "sca-ex-ch8"   "qstat"       
# [6] "m-bnd"        "m-gs1n3-5301" "sca-ex8-6"    "w-tb3n6ms"    "sp5may"      

ch08. <- paste(TsayDir, ch08[, "file"], sep="")
(ch08.datNames <- make.names(ch08[, 1]))

##
## 1.  m-ibmsp2699.txt 
##     Monthly returns of IBM and S&P 500
##

readLines(ch08.[1], 4)
m.ibmsp2699 <- read.zoo(ch08.[1], format="%Y%m%d",
      col.names=c("Date", "IBM", "SP500"))
str(m.ibmsp2699)

##
## 2.  m-ibmspln.txt
##     log returns ... 
##
readLines(ch08.[2], 4)
# ... in percent 
m.ibmspln. <- read.table(ch08.[2],
       col.names=c("IBM.log.rtn.pct", "SP500.log.rtn.pct"))
str(m.ibmspln.)

##
## 3.  m-ibmsplin.dat
##
readLines(ch08.[3], 4)
tst <- read.table(ch08.[3],
       col.names=c("IBM.log.rtn.pct", "SP500.log.rtn.pct"))
all.equal(m.ibmspln., tst)
# TRUE

# Combine m.ibmsp2699 and m.imbsplin into on zoo object
m.ibm. <- cbind(coredata(m.ibmsp2699), m.ibmspln.)
str(m.ibm.)
m.ibmTime <- index(m.ibmsp2699)
str(m.ibmTime)

m.ibmsp2699ln <- zoo(m.ibm., as.yearmon2(m.ibmTime))
str(m.ibmsp2699ln)                         
m.ibmsp2699ln[1:2, ]
index(m.ibmsp2699ln)[1:2]
                        
##
## 4.  sca-ex-ch8.txt 
##   The SCA commands used to analyze the series
##
readLines(ch08.[4], 4)

##
## 5.  qstat.f 
##     Source code of a Fortran program for multivariate Q-stat
##
readLines(ch08.[5], 4)

##
## 6.  m-bnd
##     Monthly simple returns of bond indexes
##     January 1942 - December 1999
##
readLines(ch08.[6], 4)
m.bnd. <- read.table(ch08.[6],
  col.names=c("mature30year", "mature20year",
    "mature10year", "mature5year", "mature1year"))
sapply(m.bnd., sd)
m.bnd <- zooreg(m.bnd., 1942, freq=12)
str(m.bnd)

##
## 7.  m-gs1n3-5301.txt  
##     Monthly U.S. interest rates of Example 8.6, 
##     p. 373ff
##
readLines(ch08.[7], 4)
m.gs1n3.5301. <- read.table(ch08.[7],
                col.names=c("Treasury1year", "treasury3year", "year.mo"))
m.gs1n3Date <- readLines(ch08.[7])
m.gs1n3Date[1:11]
m.gs1n3date <- substring(m.gs1n3Date, 14)
m.gs1n3date[1:11]
m.gs1n3yrmo <- as.yearmon(m.gs1n3date, "%Y.%m")
m.gs1n3yrmo[1:11]

m.gs1n3.5301 <- zoo(m.gs1n3.5301.[1:2], m.gs1n3yrmo)
m.gs1n3.5301[1:11]
str(m.gs1n3.5301)

##
## 8.  sca-ex8-6.txt
##     SCA commands used: sca-ex8-6.txt
##
readLines(ch08.[8], 4)

##
## 9.  w-tb3n6ms
##     Weekly U.S. interest rates (3-m & 6-m)
##
readLines(ch08.[9], 4)

w.tb3n6ms.t <- readLines(ch08.[9])
w.tb3n6ms.t[1:4]
w.tb3n6ms.t. <- substring(w.tb3n6ms.t, 14)
w.tb3n6ms.t.[1:4]
w.tb3n6ms.date <- as.Date(w.tb3n6ms.t., "%Y. %m. %d")
w.tb3n6ms.date[1:4]
w.tb3n6ms. <- read.table(ch08.[9],
      col.names=c("w.tb3", "w.tb6", "year", "month", "day"))
w.tb3n6ms <- zoo(w.tb3n6ms.[1:2], w.tb3n6ms.date)
w.tb3n6ms[1:4,]
str(w.tb3n6ms)

##
## 10.  sp5may
##      Log prices of SP500 index futures and shares
##      May 1993 
readLines(ch08.[10], 4)

sp5may. <- read.table(ch08.[10],
     col.names=c("logFuture", "logPrice", "dailyAvgSomething") )
dim(sp5may.)
# 7061  3
op <- par(mfrow=c(3,1))
matplot(1:7061, sp5may.[1:2], type="l")
plot(1:7061, sp5may.[[3]], type='l')
par(op)
range(sp5may.[[3]])

plot(diff(sp5may.$logFuture), type='l')
plot(diff(sp5may.$logPrice), type='l')
plot(sp5may.$dailyAvgSomething, type='l')
table(sp5may.$dailyAvgSomething)
#-0.16501 -0.14588 -0.10376 -0.08308 -0.04439 -0.04251 -0.03686  -0.0337 
#     380      379      379      375      379      332      356      379 
#-0.01418  -0.0105        0 0.000742 0.004741 0.005841 0.006489 0.012374 
#     379      347      368      379      379      355      379      379 
#0.013812 0.017094 0.018831 
#     379      379      379 
length(table(sp5may.$dailyAvgSomething))
# 19 
table(table(sp5may.$dailyAvgSomething))

May1993 <- seq(as.Date("1993/05/01"), by=1, length=31)
(MayWkEnd <- (format(May1993, "%w") %in% c(0, 6)))

library(fCalendar)
(Hol1993 <- holidayNYSE(1993))
class(Hol1993)
(MayHol <- May1993 %in% as.Date(Hol1993))
May1993. <- May1993[!(MayWkEnd | MayHol)]
# 20 days 
# Tentative conclusion:
# This may be the first differences of daily averages of something

day <- 1+c(0, cumsum(diff(sp5may.$dailyAvgSomething)!=0))
plot(day)

avgLogFuture <- tapply(sp5may.$logFuture, day, mean)
avgLogPrice <- tapply(sp5may.$logPrice, day, mean)
avgSomething <- tapply(sp5may.$dailyAvgSomething, day, mean)

plot(avgSomething[-1], diff(avgLogFuture))
plot(avgSomething[-1], diff(avgLogPrice))
# Nope.

sp5may <- cbind(sp5may., day=day)
str(sp5may)

##
## 11.  Write the data files
##

save(m.ibmsp2699ln, file="m.ibmsp2699ln.rda")

ch08.rda <- paste(ch08.datNames, "rda", sep=".")

sel8 <- c(6:7, 9:10)

for(i in sel8)
  save(list=ch08.datNames[i], file=ch08.rda[i])

save(m.ibmsp2699ln, file="m.ibmsp2699ln.rda")

############################
            
