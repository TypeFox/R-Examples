#### Create the data objects used in chapter 12 
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

# 0.1.  Tsayfiles$ch12$text 


TsayFiles$ch12$text[1:2,]
ch12 <- with(TsayFiles$ch12, text[text[, 4]=="TRUE", ])
sort(table(ch12[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch12[, "data"])
#[1] "w-gs3n1c"      "w-gs3c"        "m-sp6299"      "m-ibmspln6299"
#[5] "m-sp5-6204"    "m-geln"       

ch12. <- paste(TsayDir, ch12[, "file"], sep="")

##
## 1.  w-gs3n1c
##     Change series of weekly US interest rates (3-y & 1-y)
##     from Jan. 5, 1962, to Sep. 10, 1999 
##
readLines(ch12.[1], 4)
w.gs3n1c. <- read.table(ch12.[1],
            col.names=c("gs3", "gs1"))
str(w.gs3n1c.)
# 1966 obs. of 2 vars 
w.gs3n1Date <- seq(as.Date("19620105", "%Y%m%d"),
                   as.Date("19990919", "%Y%m%d"), 7)
str(w.gs3n1Date)
# 1968 observations
# extract from w.gs1n36299 and compare ...

data(w.gs1n36299)
w.gs1n3. <- diff(window(w.gs1n36299, start=as.Date("1962-01-05"),
    end=as.Date("1999-09-10")))[, 2:1]

all.equal(coredata(w.gs1n3.[, "gs1"]), w.gs3n1c.[, "gs1"])
# TRUE
all.equal(coredata(w.gs1n3.[, "gs3"]), w.gs3n1c.[, "gs3"])

w.gs3n1c <- w.gs1n3.[, 2:1]
str(w.gs3n1c)
w.gs3n1c[1:4,]

##
## 2.  w-gs3c
##     Change series of weekly US 3-yr interest rate
##     March 18, 1988 to Sept. 10, 1999
readLines(ch12.[2], 4)
w.gs3c. <- scan(ch12.[2])
w.gs3c <- window(w.gs3n1c[, 1], start=as.Date("1988-03-18"),
                 end = as.Date("1999-09-10"))
all.equal(coredata(w.gs3c), w.gs3c.)
# TRUE 

##
## 3.  m-sp6299
##     Monthly log returns of S&P 500 index
##
readLines(ch12.[3], 4)
m.sp6299. <- scan(ch12.[3])
# 456 items
data(m.ibmspln)
str(m.ibmspln)
m.ibmspDate <- index(m.ibmspln)
str(m.ibmspDate)
m.ibmspDate[1:2]

m.sp6299 <- window(m.ibmspln[, 2], start=yearmon(1962),
                   end=yearmon(1999+11/12))
str(m.sp6299)
# 456 items
all.equal(m.sp6299., coredata(m.sp6299))
#"Mean relative  difference: 7.446543e-05"
plot(m.sp6299.-m.sp6299)
range(m.sp6299.-m.sp6299)
m.sp6299.[1:4]
m.sp6299[1:4]
# m.sp6299 is 6 significant digits;  m.sp6299. is onlly 4

##
## 4.  m.ibmspln6299
##     Monthly log returns of IBM stock & SP 500
##
readLines(ch12.[4], 4)
m.ibmspln6299. <- read.table(ch12.[4])
# 912 items read 
m.ibmspln6299 <- window(m.ibmspln, start=yearmon(1962),
                   end=yearmon(1999+11/12))
all.equal(as.numeric(unlist(m.ibmspln6299.)),
          as.numeric(m.ibmspln6299))


##
## 5.  m.sp5.6204
##     Monthly log prices of S&P 500 index
##
readLines(ch12.[5], 4)
m.sp5.6204. <- scan(ch12.[5])
# 515 items
str(m.ibmspln)
# date range:  1926 through Dec. 1999 ... NOT Nov. 2004

m.sp6299Date <- index(m.sp6299)
range(m.sp6299Date)
table(diff(m.sp6299Date))
# OK

m.sp5.6204 <- zooreg(m.sp5.6204., 1962, freq=12)
str(m.sp5.6204)

##
## 6.  m.geln
##     Monthly log returns of GE stock 
##
readLines(ch12.[6], 4)
m.geln. <- scan(ch12.[6])
m.geln <- zooreg(m.geln., start=1926, freq=12)
str(m.geln)

##
## 7.  Write the data file
##

ch12.datNames <- make.names(ch12[, 1])
ch12.rda <- paste(ch12.datNames, "rda", sep=".")

nObj1 <- length(ch12.datNames)
for(i in 1:nObj1)
  save(list=ch12.datNames[i],
       file=ch12.rda[i])
