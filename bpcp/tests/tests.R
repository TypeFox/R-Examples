
## The following gave an error for versions<1.2.7
## there was a problem if all the censoring happended 
## at times where there was also failures
library(bpcp)
x<-c(1,1,2,3)
status<-c(0,1,1,1)
plot(c(0,5),c(0,1),type="n")
## also the default color for the survival line in 
## lines.kmciLR was white...gray(1), 
## it should be black...gray(0) 
lines(bpcp(x,status),linetype="surv")
