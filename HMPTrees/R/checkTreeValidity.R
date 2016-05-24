checkTreeValidity <-
function(data, samples=1, epsilon=0.0001, split="."){
if(missing(data))
stop("A valid data set is required.")

overAllValid <- NULL
if(samples[1] == 0)
samples <- 1:ncol(data)

for(sample in samples){
if(is.na(as.numeric(sample)))
stop(sprintf("%s must be a number", as.character(sample)))
if(sample > ncol(data))
stop(sprintf("%s is larger than the bounds of the data set", as.character(sample)))
if(sample < 1)
stop(sprintf("%s is smaller than the bounds of the data set", as.character(sample)))

tempData <- data[, sample, drop=FALSE]
nameSplit <- strsplit(rownames(tempData), split, fixed=TRUE)

for(i in 1:nrow(tempData)){ 
if(length(nameSplit[[i]]) == 1)
valid <- checkTreeValidHelp(tempData, i, epsilon, split)
}
overAllValid <- c(overAllValid, valid)
}

return(overAllValid)
}
