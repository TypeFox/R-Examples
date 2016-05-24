subsetData <-
function(data, site, region){
if(!missing(site)){
data <- subset(data, data[,1] == site)
data <- data[,-1, drop=FALSE]
}

if(!missing(region)){
data <- subset(data, data[,1] == region) 
data <- data[,-1, drop=FALSE]
}

return(data)
}
