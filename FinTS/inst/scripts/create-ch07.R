#### Create the data objects used in chapter 7 
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

# 0.1.  Tsayfiles$ch07$text 

TsayFiles$ch07$text[1:2,]
ch07 <- with(TsayFiles$ch07, text[text[, 4]=="TRUE", ])
sort(table(ch07[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch07[, "data"])

#[1] "d-ibm6298"   "example7-3a" "example7-3b" "d-intc7297"  "d-ibmln98wm"
#[6] "d-ibml25x"  

ch07. <- paste(TsayDir, ch07[, "file"], sep="")

##
## 1.  d-ibm6298
##     Daily returns in percentages of IBM stock (9190 obs)
##
readLines(ch07.[1], 4)
d.ibm6298 <- read.zoo(ch07.[1], "%Y%m%d",
                      col.names=c("Date", "IBM"))
str(d.ibm6298)

##
## 2, 3.  example7-3a/b
##     RATS programs for AR(2)-GARCH(1,1) 
##

##
## 4.  d-intc7297
##     Daily log returns of Intel stock (Example 7.4)
##
readLines(ch07.[4], 4)
d.intc7297 <- read.zoo(ch07.[4], format="%Y%m%d",
      col.names=c("Date", "Intel") )
str(d.intc7297)

##
## 5.  d-ibmln98wm
##     Mean-corrected daily log returns of IBM
##
readLines(ch07.[5], 4)
d.ibmln98wm <- read.table(ch07.[5], col.names=c("day", "meanCorrectedLogRtns") )
sapply(d.ibmln98wm, range)
#      day      IBM
#[1,]    1 -0.26132
#[2,] 9190  0.12127
table(diff(d.ibmln98wm$day))
# all 1

##
## 6.  d-ibml25x
##     The explanatory variables on page 294
##
readLines(ch07.[6], 4)
d.ibml25x <- read.table(ch07.[6], 
    col.names=c("Q4", "drop2.5pct", "nOfLast5outside2.5pct", "annualTrend",
      "GARCH1.1volatility") )
str(d.ibml25x)

# *** COMBINE d.ibm6298, d.ibmln98wm, and d.ibml25x
# *** all three have 9190 observations relating to IBM, 
# *** and apparently spanning the same time period

d.ibm6298c <- cbind(dailySimpleRtns=as.numeric(d.ibm6298),
                    d.ibmln98wm, d.ibml25x)
str(d.ibm6298c)
d.ibm6298wmx <- zoo(d.ibm6298c, index(d.ibm6298))
str(d.ibm6298wmx)

# Check Q4 = [Oct, Nov, Dec]
POSIX6298 <- as.POSIXlt(index(d.ibm6298wmx))
str(POSIX6298)
utils:::str.default(POSIX6298)
table(POSIX6298$mon, d.ibm6298wmx[,'Q4'])
# YES!  

# Check drop2.5pct
d.ibm6298wmx[1:11,]
drop2.5pct <- as.numeric(d.ibm6298wmx[-9190, "meanCorrectedLogRtns"]<=(-0.025))
str(drop2.5pct)

tst <- as.numeric(d.ibm6298wmx[-1, "drop2.5pct"])
str(tst)
table(drop2.5pct, tst)
# checks 

# check nOfLast5outside2.5pct

str(d.ibm6298)


n5s <- n5 <- rep(0, 9185)
seq0 <- 6:9190
sDR2.5 <- as.numeric(abs(d.ibm6298wmx[, "dailySimpleRtns"])>=0.025)
mclr2.5 <- as.numeric(abs(d.ibm6298wmx[, "meanCorrectedLogRtns"])>=0.025)
all.equal(sDR2.5, mclr2.5)

for(i in 1:5){
  seqi <- seq0-i
  n5s <- (n5s+sDR2.5[seqi])
  n5 <- (n5+mclr2.5[seqi])
}

tst <- as.numeric(d.ibm6298wmx[-(1:5), "nOfLast5outside2.5pct"])
table(n5s, tst)
# close but ... 
table(n5, tst)
# TRUE

# check annual trend
table(as.POSIXlt(index(d.ibm6298wmx))$year)
table(aT)
aT <- (as.POSIXlt(index(d.ibm6298wmx))$year-61)/38
table(aT)

plot(aT, as.numeric(d.ibm6298wmx[, "annualTrend"]))
sum(diag(table(aT, as.numeric(d.ibm6298wmx[, "annualTrend"]))))
# 9190 Good  
all.equal(aT, as.numeric(d.ibm6298wmx[, "annualTrend"]))
# mean rel diff = 5e-5 round off         

##
## 7.  Write the data file
##

ch07.datNames <- make.names(ch07[, 1])
ch07.rda <- paste(ch07.datNames, "rda", sep=".")

#sel7 <- c(1, 4:6)
#for(i in sel7)save(list=ch07.datNames[i], file=ch07.rda[i])

save(d.ibm6298wmx, file="d.ibm6298wmx.rda")
save(d.intc7297, file="d.intc7297.rda")
