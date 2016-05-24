#### Create the data objects used in chapter 9 
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

# 0.1.  Tsayfiles$ch09$text 

TsayFiles$ch09$text[1:2,]

ch09 <- with(TsayFiles$ch09, text[text[, 4]=="TRUE", ])
sort(table(ch09[, "data"]))
# Confirm:  all unique
# Exercises or text in other chapters may use
# some of these same data;  we need to check then.  

as.vector(ch09[, "data"])
#[1] "m-fac9003"        "m-cpice16-dp7503" "m-barra-9003"     "m-5cln"          
#[5] "m-bnd"            "m-apca0103"
                                        
ch09. <- paste(TsayDir, ch09[, "file"], sep="")
(ch09.datNames <- make.names(ch09[, 1]))

##
## 1.  m-fac9003 
##     Monthly stock returns of Table 9.1
##
readLines(ch09.[1], 4)
m.fac9003. <- read.table(ch09.[1], 
   col.names=c("AA", "AGE", "CAT", "F", "FDX", "GM", "HPQ",
     "KMB", "MEL", "NYT", "PG", "TRB", "TXN", "SP5"))
str(m.fac9003.)
(m.facYrMo <- yearmon(1990+(0:167)/12))
round(sapply(m.fac9003., mean), 2)
round(sapply(m.fac9003., sd), 2)

m.fac9003 <- zoo(m.fac9003., m.facYrMo)
str(m.fac9003)
##
## 2.  m-cpice16-dp7503
##     Monthly macroeconomic variables: (CPI & CE16;  p. 412)
##
readLines(ch09.[2], 4)
m.cpice16.dp7503. <- read.table(ch09.[2], col.names=c("CPI", "CE16"))
str(m.cpice16.dp7503.)
m.cpiYrMo <- yearmon(1975+(0:335)/12)
m.cpice16.dp7503 <- zoo(m.cpice16.dp7503., m.cpiYrMo)
str(m.cpice16.dp7503)
##
## 3.  m-barra-9003 
##     Monthly excess returns of Table 9.2 (p. 416) 
##
readLines(ch09.[3], 4)
m.barra.9003. <- read.table(ch09.[3], col.names=c(
     "AGE", "C", "MWD", "MER", "DELL", "HPQ", "IBM", "AA", "CAT", "PG"))
str(m.barra.9003.)
m.barraYrMo <- yearmon(1990+(0:167)/12)

m.barra.9003 <- zoo(m.barra.9003., m.barraYrMo)

##
## 4.  m-5cln 
##     Monthly log returns, in percentages, of
##     IBM, HPQ, INTC, MER & MWD stocks
##
readLines(ch09.[4], 4)
m.5cln. <- read.table(ch09.[4], col.names=c(
       "IBM", "HPQ", "INTC", "MER", "MWD"))
str(m.5cln.)
m.5clnYrMo <- yearmon(1990+(0:119)/12)
m.5cln <- zoo(m.5cln., m.5clnYrMo)

##
## 5.  m-bnd
##     Monthly returns of U.S. bond indices
##
readLines(ch09.[5], 4)
m.bnd[1:4,]
tst <- read.table(ch09.[5])
for(i in 1:5)
  cat(all.equal(tst[[i]], coredata(m.bnd[, i])), "")
#TRUE TRUE TRUE TRUE TRUE 
##
## 6.  m-apca0103.txt
##     Monthly returns of 40 stocks in Table 9.6:
##     (Company ID, date, return)
##
readLines(ch09.[6], 4)
m.apca0103. <- readLines(ch09.[6])
CoID <- substring(m.apca0103., 1, 7)
CoID[1:4]
m.apcaDate <- as.Date(substring(m.apca0103., 8, 20), "%Y%m%d")
m.apcaDate[1:4]
(oops <- which(diff(m.apcaDate)<=0))
m.apca0103.[35:37]
table(m.apcaDate)
plot(m.apcaDate)
m.apca0103.[c(1, 36)]

table(format(m.apcaDate[1:36], "%w"))

m.apca0103 <- data.frame(CompanyID=as.numeric(CoID),
                         date=m.apcaDate,
                         return=as.numeric(substring(m.apca0103., 21)))
str(m.apca0103)
table(m.apca0103$CompanyID)

##
## 7.  Write the data files
##

ch09.rda <- paste(ch09.datNames, "rda", sep=".")

sel9 <- c(1:4, 6)

for(i in sel9)
  save(list=ch09.datNames[i], file=ch09.rda[i])
