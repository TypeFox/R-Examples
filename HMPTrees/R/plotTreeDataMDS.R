plotTreeDataMDS <-
function(data, main="Tree MDS Comparisons", calcMLE=TRUE, mleTitles, 
dotColors=c("red", "orange", "blue", "green", "yellow", "purple"), dotSizes=c(1, 2), 
showNames=FALSE, returnCoords=FALSE){
normalDotSizes <- NULL
mleDotSizes <- NULL
normalColors <- NULL
mleColors <- NULL
mleLabels <- NULL

if(missing(data))
stop("At least 1 valid data set is required.")

if(class(data) == "data.frame" || class(data) == "matrix"){ #turn a single dataset into a list
tempData <- data
data <- NULL
data[[1]] <- tempData
}

if(length(dotColors) < length(data))
dotColors <- rainbow(length(data))

twoColors <- length(dotColors) >= length(data)*2 #2 colors per data set
dataCount <- length(data)
titles <- !missing(mleTitles)
if(titles) #make sure we have the same number of titles as data sets
titles <- length(mleTitles) == length(data)

for(i in 1:dataCount){
tempData <- data[[i]]
if(twoColors){
colorLoc <- i*2-1
colorLoc2 <- colorLoc+1
}else{
colorLoc <- i
colorLoc2 <- i
}

normalColors <- c(normalColors, rep(dotColors[colorLoc], ncol(tempData)))
mleColors <- c(mleColors, dotColors[colorLoc2])
normalDotSizes <- c(normalDotSizes, rep(dotSizes[1], ncol(tempData)))
mleDotSizes <- c(mleDotSizes, dotSizes[2])
if(calcMLE && titles)
mleLabels <- c(mleLabels, mleTitles[i])
}

mData <- mergeDataSets(data, calcMLE)
tData <- t(mData)
loc <- cmdscale(dist(tData), k=2)
x <- loc[,1]
y <- -loc[,2]

colors <- c(normalColors, mleColors, "black")
sizes <- c(normalDotSizes, mleDotSizes, dotSizes[2])
plot(x, y, pch=19, xlab="MDS 1", ylab="MDS 2", pty="s", col=colors, cex=sizes, main=main)

if(calcMLE && !showNames){ #only plot mle titles if we are calculating them
if(dataCount > 1){ #dont add a combined mle if we dont make one
if(titles){
mleTitles <- c(mleTitles, "Combined MLE")
}else{
mleTitles <- paste("Data", 1:length(data), "MLE")
}
dataCount <- dataCount + 1
}
text(x[(nrow(tData)-dataCount+1):nrow(tData)], y[(nrow(tData)-dataCount+1):nrow(tData)], mleTitles, pos=3, cex=.75)
}

if(showNames) #Plot the names of the samples
text(x, y, colnames(mData), pos=3, cex=.75)

if(returnCoords)
return(list(x=x, y=y))
}
