uncertainty <- function(index, tocd){ 

thist <- hist(unique(index), breaks=c(sort(tocd$Threshold)), plot=FALSE)
tocd$counts <- c(0, thist$counts[length(thist$counts):1], 0)

# calculate uncertainty
uncertain <- 0
for (i in 2:(nrow(tocd)-2)){
  if(tocd[i,"counts"] > 1) {
    area <- (tocd[i,"falseAlarms1"] - tocd[i-1,"falseAlarms1"])*(tocd[i,"Model1"] - tocd[i-1,"Model1"])
    uncertain <- uncertain + area
  }
}

i <- i+1
if(tocd[i,"counts"] > 2) {
  area <- (tocd[i,"falseAlarms1"] - tocd[i-1,"falseAlarms1"])*(tocd[i,"Model1"] - tocd[i-1,"Model1"])
  uncertain <- uncertain + area
}

return(uncertain)
}