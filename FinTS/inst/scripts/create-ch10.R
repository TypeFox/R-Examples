#### Create the data objects used in chapter 10 
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

# 0.1.  Tsayfiles$ch10$text 

TsayFiles$ch10$text[1:2,]

ch10 <- with(TsayFiles$ch10, text[text[, 4]=="TRUE", ])
sort(table(ch10[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch10[, "data"])

# [1] "d-hkja"        "hkja-c"        "hkja-c1"       "m-pfe6503"    
# [5] "m-mrk6503"     "m-ibmsp2699"   "ibmsp-ex92"    "ibmsp-ex92q"  
# [9] "ibmsp-choles"  "d-spcscointc"  "cholesky-ex93"

ch10. <- paste(TsayDir, ch10[, "file"], sep="")
(ch10.datNames <- make.names(ch10[, 1]))

##
## 1.  d.hkja
##     Daily log returns of HK and Japan market indices
##     (Example 10.1):  Data file (491 data pts)
##

readLines(ch10.[1], 4)
d.hkja. <- read.table(ch10.[1],
                      col.names=c("HongKong", "Japan"))
str(d.hkja.)
plot(d.hkja.$HongKong, type="l")
plot(d.hkja.$HongKong[1:469], type="l")

d.hkjaDate <- seq(as.Date("19960101", "%Y%m%d"),
                  by=1, length=491)

d.hkja <- zoo(d.hkja., d.hkjaDate)
str(d.hkja)
##
## 2-3.  hkja-c.rats, hkja-c1.rats
##     RATS programs 
##

##
## 4.  m-pfe6503
##     Monthly simple returns of Pfizer stocks
##     (including dividends) 
##
readLines(ch10.[4], 4)
m.pfe6503 <- read.yearmon(ch10.[4], format="%Y%m%d",
      col.names=c("date", "return"))

##
## 5.  m-mrk6503
##     Monthly simple returns of Merck stocks
##     (including dividends)
##
readLines(ch10.[5], 4)
m.mrk6503 <- read.yearmon(ch10.[5], format="%Y%m%d",
       col.names=c("date","return"))

##
## 6.  m-ibmsp2699
##     Monthly returns of IBM and S&P 500
##
readLines(ch10.[6], 4)
m.ibmsp2699 <- read.yearmon(ch10.[6], format="%Y%m%d",
           col.names=c("date", "IBM", "SP"))
m.ibmsp2699[1:3,]
m.ibmsp2699ln[1:3,]
all.equal(m.ibmsp2699, m.ibmsp2699ln[, 1:2])
# problems ...
str(m.ibmsp2699)
str(m.ibmsp2699ln)
all.equal(coredata(m.ibmsp2699), coredata(m.ibmsp2699ln[, 1:2]))
all.equal(as.numeric(m.ibmsp2699[,1]),
          as.numeric(m.ibmsp2699ln[, 1]))
# TRUE
all.equal(as.numeric(m.ibmsp2699[,2]),
          as.numeric(m.ibmsp2699ln[, 2]))
# TRUE
# OK

##
## 7-9.  ibmsp-ex92.rats, ibmsp-ex92q.rats, ibmsp-choles.rats
##     RATS code for GARCH & Choleski
##

##
## 10.  d-spcscointc.txt
##      Daily log returns of S&P 500, Cisco and Intel stocks
##      Data (3 columns)
##
readLines(ch10.[10], 4)
d.spcscointc <- read.table(ch10.[10], col.names=
             c("SP500", "Cisco", "Intel"))
str(d.spcscointc)
#  2275 obs. of 3 vars
plot(d.spcscointc$SP500, type="l")
plot(d.spcscointc$Cisco, type="l")
plot(d.spcscointc$Intel, type="l")

d.spcscoDates <- seq(as.Date("19910102", "%Y%m%d"),
                     as.Date("19991231", "%Y%m%d"), 1) 
str(d.spcscoDates)
# 3286 obs.
d.spcscoWkEnd <- (format(d.spcscoDates, "%w") %in% c(0, 6))

library(fCalendar)
d.spcscoHol <- (d.spcscoDates %in%
                as.Date(holidayNYSE(1991:1999)))
d.spcscoDts <- d.spcscoDates[!(d.spcscoWkEnd | d.spcscoHol)]
str(d.spcscoDts)
# 2280
# This is different from the number of observations
# in the series, 2275;  leave it without dates for now.  

##
## 11.  Write data files
##
ch10.rda <- paste(ch10.datNames, "rda", sep=".")

sel10 <- c(1, 4:5, 10)

for(i in sel10)
  save(list=ch10.datNames[i], file=ch10.rda[i])

#######################
# fix changes made after writing as.yearmon2 and read.yearmon

m.obj10 <- 4:5
for(i in m.obj10)
  save(list=ch10.datNames[i], file=ch10.rda[i])
