#### Create the data objects used in chapter 2 
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

# 0.1.  Tsayfiles$ch02$text 

TsayFiles$ch02$text[1:2,]
ch02 <- with(TsayFiles$ch02, text[text[, 4]=="TRUE", ])
sort(table(ch02[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch02[, "data"])

# [1] "m-ibm2697"    "m-vw2697"     "q-gnp4791"    "m-ibm3dx2603" "m-3m4697"    
# [6] "q-gdp4703"    "d-sp9003lev"  "q-jnj"        "m-decile1510" "w-gs1n36299" 

ch02. <- paste(TsayDir, ch02[, "file"], sep="")

##
## 1.  m-ibm2697
##     Monthly IBM stock returns
##

readLines(ch02.[1], 4)
m.ibm2697a <- scan(ch02.[1])
# read 864 items
# (97-25)*12 = 864
m.ibm2697 <- zooreg(m.ibm2697a, 1926, freq=12)
m.ibm2697[1:22]

##
## 2.  m-vw2697
##     Monthly returns of VW index
##
readLines(ch02.[2], 4)
m.vw2697a <- scan(ch02.[2])
# read 864 items
# (97-25)*12 = 864
m.vw2697 <- zooreg(m.vw2697a, 1926, freq=12)
m.vw2697[1:22]

##
## 3.  q.gnp4791
##     Growth rate of U.S. quarterly real gnp
##
readLines(ch02.[3], 4)
q.gnp4791a <- scan(ch02.[3])
# read 176 items
# (91-47)*4 = 864
q.gnp4791 <- zooreg(q.gnp4791a, c(1947, 2), freq=4)
q.gnp4791[1:22]

##
## 4.  m-ibm3dx2603
##     Monthly returns of EW index: (Date, IBM, VW, EW & SP)
##
readLines(ch02.[4], 4)
m.ibm3dx2603 <- read.yearmon(ch02.[4], format="%Y%m%d",
      col.names=c("Date", "IBM", "VW", "EW", "SP") )
m.ibm3dx2603[1:4,]

##
## 5.  m.3m4697
##     Monthly simple returns of 3M stock
##
readLines(ch02.[5], 4)
m.3m4697 <- read.yearmon(ch02.[5], format="%Y%m%d",
      col.names=c("Date", "MMM") )
m.3m4697[1:4,]

##
## 6.  q.gdp4703
##     U.S. quarterly GDP
##
readLines(ch02.[6], 4)
q.gdp4703a <- readLines(ch02.[6])
q.gdp4703a[1:2]
q.gdp4703d <- as.yearmon(q.gdp4703a, "%Y %m")
q.gdp4703d[1:4]
class(q.gdp4703d)

q.gdp4703q <- as.yearqtr(q.gdp4703d)
q.gdp4703q[1:2]
names(q.gdp4703q) <- q.gdp4703d
q.gdp4703q[1:2]
index(q.gdp4703q)[1:2]

q.gdp4703 <- zoo(as.numeric(substring(q.gdp4703a, 8)),
                 q.gdp4703q)
q.gdp4703[1:4]
#q.gdp4703t <- strptime(q.gdp4703a, "%Y %m")
#q.gdp4703t[1:4]
# NA 

##
## 7.  d.sp9003lev
##     Daily values of S&P 500 index
##
readLines(ch02.[7], 4)
d.sp9003lev <- read.zoo(ch02.[7], format="%Y%m%d",
                        col.names=c("date", "sp"))
d.sp9003lev[1:4]

##
## 8.  q.jnj
##     Quarterly earnings of JNJ (1960-1980)
readLines(ch02.[8], 4)
q.jnj0 <- scan(ch02.[8])
# read 84 items
# (80-59)*4
q.jnj0[1:4]
q.jnj <- zooreg(q.jnj0, 1960, freq=4)
q.jnj[1:4]

##
## 9.  m.decile1510
##     Monthly simple returns of Deciles 1, 5, 10: m-decile1510.txt
##
readLines(ch02.[9], 4)
m.decile1510 <- read.yearmon(ch02.[9], format="%Y%m%d",
     col.names=c("date", "Decile1", "Decile5", "Decile10") )       
m.decile1510[1:4,]

##
## 10.  w.gs1n36299
##      Weekly 1-yr & 3-yr interest rates
##
readLines(ch02.[10], 4)
#w.gs1n36299 <- read.zoo(ch02.[10], ...) ERRORS IN DATES;

# *** USE CORRRECTED FILE as follows:  
w.gs1n36299File <- paste(TsayDir, "w-gs1n36299r.txt", sep="")
w.gs1n36299 <- read.zoo(w.gs1n36299File, format="%Y-%m-%d", header=TRUE)
w.gs1n36299[1:4,]

##
## 11.  Write the data file
##

ch02.datNames <- make.names(ch02[, 1])
ch02.rda <- paste(ch02.datNames, "rda", sep=".")

nObj2 <- length(ch02.datNames)
for(i in 1:nObj2)
  save(list=ch02.datNames[i], file=ch02.rda[i])

##########################
m.obj2 <- c(4:5, 9)
for(i in m.obj2)
  save(list=ch02.datNames[i], file=ch02.rda[i])
##############################
  save(list=ch02.datNames[6], file=ch02.rda[6])
