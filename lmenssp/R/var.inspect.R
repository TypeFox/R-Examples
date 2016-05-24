var.inspect <- function(resid, timeVar, binwidth, numElems = 0, irregular = T){

if(irregular == TRUE){

data     <- data.frame(time = timeVar, resid = resid)
var.base <- var(data[data$time == 0, "resid"])
num.base <- length(data[data$time == 0, "resid"])
bin.lims <- seq(0.00000001, round(max(data$time)), by = binwidth)

var.save <- nums <- c()

for(i in 2 : length(bin.lims)){
var.save <- c(var.save, var(data[data[, "time"] >= bin.lims[i - 1] & data[, "time"] < bin.lims[i], "resid"]))
nums     <- c(nums, length(data[data[, "time"] >= bin.lims[i - 1] & data[, "time"] < bin.lims[i], "resid"]))
}

bin.mids <- c()
for (i in 2 : length(bin.lims)){
bin.mids <- c(bin.mids, (bin.lims[i-1] + bin.lims[i])/2)
}

bin.mids <- c(0, bin.mids)
var.save <- c(var.base, var.save)
nums     <- c(num.base, nums)

plot(bin.mids[which(nums > numElems)], 
     var.save[which(nums > numElems)], 
     ylim = c(min(var.save[which(nums > numElems)], na.rm = T) - 0.05, max(var.save[which(nums > numElems)], na.rm = T) + 0.05),
     type = "l", ylab = "Variances", xlab = "Time")

output <- list()
output$bin.mids   <- bin.mids
output$variances  <- var.save
output$bin.sizes  <- nums

}

if(irregular == FALSE){

time.save <- sort(unique(timeVar))
var.save  <- tapply(resid, timeVar, function(x) var(x))
size.save <- tapply(resid, timeVar, function(x) length(x))

plot(time.save, var.save, ylim = c(min(var.save) - 0.1, max(var.save) + 0.1),
     type = "l", ylab = "Variances", xlab = "Time")

output           <- list()
output$time      <- time.save
output$variances <- var.save
output$sizes     <- size.save

}

output

}
