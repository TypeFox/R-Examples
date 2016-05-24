readDates3to1 <- function(file, YMD=c('Year', 'Month', 'Day'), 
                          ...){
##  
##  1.  read.csv
##  
  dat <- read.csv(file, ...)
##
## 2.  Dates3to1
##
  Dates3to1(dat, YMD)
}
