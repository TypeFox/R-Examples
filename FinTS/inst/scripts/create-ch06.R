#### Create the data objects used in chapter 6 
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

# 0.1.  Tsayfiles$ch06$text 

TsayFiles$ch06$text[1:2,]

ch06 <- with(TsayFiles$ch06, text[text[, 4]=="TRUE", ])
sort(table(ch06[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch06[, "data"])

#[1] "m-intc7303"      "exch-perc"       "sp500"           "m-ibm2697"      
#[5] "d-ibmvwewsp6203" "m-ibmspln"       "m-ibmsplnsu"     "d-sp8099"       

ch06. <- paste(TsayDir, ch06[, "file"], sep="")
(ch06.datNames <- make.names(ch06[, 1]))

##
## 1.  d-ibmy98
##     Daily simple returns of IBM stock in 1998
##
readLines(ch06.[1], 4)
d.ibmy98 <- read.zoo(ch06.[1], format="%Y%m%d",
                     col.names=c("Date", "IBM"))
str(d.ibmy98)

## 2.  d-cscoy99
##     Daily log returns of Cisco stock in 1999
##
readLines(ch06.[2], 4)
d.cscoy99 <- read.zoo(ch06.[2], format="%Y%m%d",
                      col.names=c("Date", "Cisco"))
str(d.cscoy99) 

##
## 3.  Write the data files
##

ch06.rda <- paste(ch06.datNames, "rda", sep=".")

for(i in 1:2)
  save(list=ch06.datNames[i], file=ch06.rda[i])
