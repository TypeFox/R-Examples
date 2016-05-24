#### Create the data objects used in chapter 4 
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

# 0.1.  Tsayfiles$ch04$text 

TsayFiles$ch04$text[1:2,]

ch04 <- with(TsayFiles$ch04, text[text[, 4]=="TRUE", ])
sort(table(ch04[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch04[, "data"])

#[1] "m-unrate"        "d-ibmvwewsp6203" "m-3m4697"        "star"           
#[5] "q-gnp4791"       "w-3mtbs7097"     "m-ibmln2699"     "q-unemrate"     
#[9] "nnet-ibm"       

ch04. <- paste(TsayDir, ch04[, "file"], sep="")
(ch04.datNames <- make.names(ch04[, 1]))

##
## 1.  m-unrate
##     Monthly U.S. civilian unemployment rate(48-04)
##

readLines(ch04.[1], 4)
m.unrate. <- readLines(ch04.[1])
m.unrate <- zoo(as.numeric(substring(m.unrate., 11)),
                as.yearmon2(as.Date(m.unrate., "%Y %m %d")))
m.unrate[1:4]
index(m.unrate)[1:4]

##
## 2.  d-ibmvwewsp6203
##     Daily returns of IBM stock
##
readLines(ch04.[2], 4)
tst.ibmvwewsp6203 <- read.zoo(ch04.[2], format="%Y%m%d",
          col.names=c("date", "IBM", "VW", "EW", "SP") )
data(d.ibmvwewsp6203)
all.equal(tst.ibmvwewsp6203, d.ibmvwewsp6203)
# TRUE

##
## 3.  m-3m4697
##     Monthly simple returns of 3M stock
##
readLines(ch04.[3], 4)
tst.3m4697 <- read.zoo(ch04.[3], format="%Y%m%d",
                       col.names=c("date", "MMM"))
data(m.3m4697)
all.equal(tst.3m4697, m.3m4697)
# TRUE 

##
## 4.  star.rats
##     RATS program for smooth TAR:  ignore
##
(star.rats <- readLines(ch04.[4]))

##
## 5.  q-gnp4791
##     Quarterly growth rates of U.S. gnp
##
readLines(ch04.[5], 4)
tst.gnp4791 <- scan(ch04.[5])

data(q.gnp4791)
str(q.gnp4791)
all.equal(tst.gnp4791, as.numeric(q.gnp4791))
# TRUE

##
## 6.  w-3mtbs7097
##     Weekly 3-month Treasury Bill rates
##
readLines(ch04.[6], 4)
w.3mtbs7097. <- readLines(ch04.[6])

w.3mtbs7097 <- zoo(as.numeric(substring(w.3mtbs7097., 11)),
                   as.Date(w.3mtbsDate., "%Y %m %d"))
w.3mtbs7097[1:4]
range(index(w.3mtbs7097))

##
## 7.  m-ibmln2699
##     Monthly log returns, in percentages, of IBM stock
##
readLines(ch04.[7], 4)
m.ibmln2699. <- scan(ch04.[7])
# read 888 items
# (99-25)*12 = 888 
m.ibmln2699 <- zooreg(m.ibmln2699., 1926, freq=12)
m.ibmln2699[1:4]

##
## 8.  q-unemrate
##     Quarterly unemployment rates
##
readLines(ch04.[8], 4)
q.unemrate. <- scan(ch04.[8] )
# read 184 items
# (1993-1947)*4 = 184

q.unemrate <- zooreg(q.unemrate., 1948, freq=4)
q.unemrate[1:4]
##
## 9.  nnet-ibm.sor
##     R and S commands for Example 4.7
##

##
## 10.  Write data files
##
ch04.rda <- paste(ch04.datNames, "rda", sep=".")

sel4 <- c(1, 6:8)

for(i in sel4)
  save(list=ch04.datNames[i], file=ch04.rda[i])

#######################
# correct the one file changed after writing as.yearmon2 and read.yearmon:
#    save(list=ch04.datNames[1], file=ch04.rda[1])
