graphNetworkPlot <-
function(data, type, main="Network Plot", labels, groupCounts, groupLabels){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

#Set up plot dot colors
if(missing(groupCounts)){
myColors <- "red"
}else{
allColors <- rainbow(length(groupCounts)+1)
myColors <- NULL
for(i in 1:length(groupCounts))
myColors <- c(myColors, rep(allColors[i], groupCounts[i]))
}

#Take only the first column of data if it is multi columned
if(class(data) == "data.frame" || class(data) == "matrix")
data <- data[,1]

y <- vec2mat(data, type)
g <- network::network(as.matrix(y), directed=FALSE)
if(!missing(labels))
network::network.vertex.names(g) <- as.data.frame(labels)

network::plot.network(g, mode="circle", vertex.col=myColors, label=network::network.vertex.names(g), main=main, edge.col="black")

if(!missing(groupLabels))
legend("topright", legend=groupLabels, fill=allColors, horiz=FALSE)
}
