#### Create the data objects used in chapter 11 
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

# 0.1.  Tsayfiles$ch11$text 


TsayFiles$ch11$text
ch11 <- with(TsayFiles$ch11, text[text[, 4]=="TRUE",, drop=FALSE])
sort(table(ch11[, "data", drop=FALSE]))
# Confirm:  all unique.  
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.
#
# With only 1 file, there can't be duplicates.
# However, this script is modified from scripts for other chapters
# with multiple files ... .

as.vector(ch11[, "data"])
#"aa-3rv"

ch11. <- paste(TsayDir, ch11[, "file"], sep="")

##
## 1.  aa-3rv
##     Daily realized volatility series of Alcoa stock: (5m, 10m, 20m)
##
readLines(ch11.[1], 4)

aa.3rv. <- read.table(ch11.[1], col.names=c("5m", "10m", "20m"))
str(aa.3rv.)
# 340 obs. on 3 vars.
aa.3Dates <- seq(as.Date("20030102", "%Y%m%d"),
                 as.Date("20040507", "%Y%m%d"), by=1)
aa3.WkEnd <- (format(aa.3Dates, "%w") %in% c(0, 6))

library(fCalendar)
aa.3hol <- (aa.3Dates %in% as.Date(holidayNYSE(2003:2004)))

aa.3Dt <- aa.3Dates[!(aa3.WkEnd | aa.3hol)]
str(aa.3Dt)
# 340 :-)
aa.3rv <- zoo(aa.3rv., aa.3Dt)
str(aa.3rv)

##
## 2.  Write the data file
##

ch11.datNames <- make.names(ch11[, 1])
ch11.rda <- paste(ch11.datNames, "rda", sep=".")
#save(list=ch11.datNames, file="ch11.rda")

nObj1 <- length(ch11.datNames)
for(i in 1:nObj1)
  save(list=ch11.datNames[i],
       file=ch11.rda[i])
