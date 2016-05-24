##############################################
# Section 7.3 Individual or Combined Estimates
##############################################

 library(LearnBayes)
 data(hearttransplants)
 attach(hearttransplants)

 plot(log(e), y/e, xlim=c(6,9.7), xlab="log(e)", ylab="y/e")
 text(log(e),y/e,labels=as.character(y),pos=4)
