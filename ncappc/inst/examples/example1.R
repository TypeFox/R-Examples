## Test example for ncappc 

library(ncappc)

# This example checks the processing part of time and concentration data

ncappc.db <- data.frame(TIME=c(0:10),CONC=c(0,1,2,3.5,3.6,3.2,2.4,1.6,1.3,1.1,0.8))

TIME <- ncappc.db$TIME
CONC <- ncappc.db$CONC
TIME.positive <- ncappc.db[ncappc.db$TIME >= 0,"TIME"]
CONC.positive <- ncappc.db[ncappc.db$CONC >= 0,"CONC"]


