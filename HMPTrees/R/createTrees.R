createTrees <-
function(data, samples=1, level="genus"){
if(missing(data))
stop("A valid data set is required.")

if(samples[1] == 0)
samples <- 1:ncol(data)

data <- data[order(rownames(data)),]

allTrees <- list()
i <- 1
for(sample in samples){
if(is.na(as.numeric(sample)))
stop(sprintf("%s must be a number", as.character(sample)))
if(sample > ncol(data)) 
stop(sprintf("%s is larger than the bounds of the data set", as.character(sample)))
if(sample < 1)
stop(sprintf("%s is smaller than the bounds of the data set", as.character(sample)))

oneSamp <- data[, sample, drop=FALSE]

if(sum(oneSamp) <= 0)#skips entries without data
next
tempTree <- traverseTree(oneSamp, level)
allTrees[[i]] <- tempTree
i <- i + 1
}

names(allTrees) <- colnames(data)[samples]

return(allTrees)
}
